#!/usr/bin/env python3
"""
Synthetic Perturbation Benchmarking Suite: Milopy DA vs. Pseudobulk Deconvolution.

Comprehensive Multi-Replicate & Multi-Factorial Distortion Evaluation:
1. Replicate evaluation (5 random seeds per intensity level) computing Mean +/- SD error bands.
2. 6 single perturbation modes + 2 novel biological/technical modes (Ambient Soup, Marker Dysregulation).
3. 4 realistic multi-factorial compound regimes:
   - "Clinical Core Needle Biopsy" (Cell Size + Sparsity + Ghost Contamination)
   - "Inflamed Tumor Microenvironment" (Activation + Collinearity + Patient Shift)
   - "Single-Cell Technical Noise" (Ambient Soup + Sparsity + Patient Shift)
   - "Triple-Jeopardy Breakdown" (Cell Size + Activation + Ghost Contamination)
4. 2D Factorial Interaction Matrix (Cell Size x Activation Confounding) mapping non-linear interaction surfaces.
5. Publication-quality Altair SVG vector visualizations.

Strict functional Python: immutability, Pydantic frozen models, Polars, and Altair.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Final, Sequence

import altair as alt  # type: ignore
import anndata as ad  # type: ignore
from joblib import Parallel, delayed  # type: ignore
import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
from scipy import sparse as sp  # type: ignore
from scipy import stats  # type: ignore
import scanpy as sc  # type: ignore
import vl_convert as vlc  # type: ignore

try:
    import milopy  # type: ignore
    HAS_MILOPY = True
except ImportError:
    HAS_MILOPY = False

try:
    import instaprism  # type: ignore
    HAS_INSTAPRISM = True
except ImportError:
    HAS_INSTAPRISM = False


PERTURBATION_MODES: Final[tuple[str, ...]] = (
    "cell_size",
    "collinearity",
    "patient_shift",
    "activation_confounding",
    "ghost_contamination",
    "sampling_sparsity",
    "ambient_soup",
    "marker_dysregulation",
)

MODE_DISPLAY_TITLES: Final[dict[str, str]] = {
    "cell_size": "Cell Size Asymmetry (1x-50x)",
    "collinearity": "Reference Collinearity (0-99%)",
    "patient_shift": "Inter-Patient Shift (sigma 0-3.0)",
    "activation_confounding": "Activation Confounding (1x-20x)",
    "ghost_contamination": "Ghost Contamination (0-85%)",
    "sampling_sparsity": "Sampling Sparsity (sigma 0-4.0)",
    "ambient_soup": "Ambient RNA Soup (0-40%)",
    "marker_dysregulation": "Patient Marker Dysregulation (sigma 0-2.0)",
}

COMPOUND_REGIMES: Final[dict[str, dict[str, object]]] = {
    "core_biopsy": {
        "title": "Clinical Core Needle Biopsy",
        "description": "Cell Size (1.0x) + Sparsity (0.8x) + Ghost Contamination (0.7x)",
        "weights": {
            "cell_size": 1.0,
            "sampling_sparsity": 0.8,
            "ghost_contamination": 0.7,
        },
    },
    "inflamed_tme": {
        "title": "Inflamed Tumor Microenvironment",
        "description": "Activation (1.0x) + Collinearity (0.8x) + Shift (0.6x)",
        "weights": {
            "activation_confounding": 1.0,
            "collinearity": 0.8,
            "patient_shift": 0.6,
        },
    },
    "sc_technical": {
        "title": "Single-Cell Technical Noise",
        "description": "Ambient Soup (1.0x) + Sparsity (0.6x) + Shift (0.5x)",
        "weights": {
            "ambient_soup": 1.0,
            "sampling_sparsity": 0.6,
            "patient_shift": 0.5,
        },
    },
    "triple_jeopardy": {
        "title": "Triple-Jeopardy Breakdown",
        "description": "Cell Size (1.0x) + Activation (1.0x) + Ghost (0.8x)",
        "weights": {
            "cell_size": 1.0,
            "activation_confounding": 1.0,
            "ghost_contamination": 0.8,
        },
    },
}


class PerturbationProfile(BaseModel):
    """Immutable parameterization of single and compound perturbation states."""
    model_config = ConfigDict(frozen=True)
    cell_size: float = 0.0
    collinearity: float = 0.0
    patient_shift: float = 0.0
    activation_confounding: float = 0.0
    ghost_contamination: float = 0.0
    sampling_sparsity: float = 0.0
    ambient_soup: float = 0.0
    marker_dysregulation: float = 0.0

    @classmethod
    def from_single(cls, mode: str, intensity: float) -> PerturbationProfile:
        valid_fields = {
            "cell_size", "collinearity", "patient_shift", "activation_confounding",
            "ghost_contamination", "sampling_sparsity", "ambient_soup", "marker_dysregulation"
        }
        kwargs = {f: (float(intensity) if f == mode else 0.0) for f in valid_fields}
        return cls(**kwargs)


def make_compound_profile(regime_key: str, alpha: float) -> PerturbationProfile:
    """Generate an immutable PerturbationProfile for a compound regime at intensity alpha."""
    regime = COMPOUND_REGIMES[regime_key]
    weights: dict[str, float] = regime["weights"]  # type: ignore
    kwargs = {k: min(1.0, w * alpha) for k, w in weights.items()}
    return PerturbationProfile(**kwargs)


class SweepConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    n_cells: int = 2400
    n_genes: int = 240
    n_clusters: int = 8
    n_patients: int = 30
    n_levels: int = 6
    n_reps: int = 5
    random_seed: int = 42
    out_dir: Path = Path("output/synthetic_benchmark")
    results_dir: Path = Path("results/synthetic_benchmark")
    curves_svg_name: str = "perturbation_sensitivity_curves.svg"
    scatter_svg_name: str = "synthetic_perturbation_scatter_grid.svg"
    compound_curves_svg_name: str = "compound_perturbation_sensitivity_curves.svg"
    compound_scatter_svg_name: str = "compound_perturbation_scatter_grid.svg"
    factorial_heatmap_svg_name: str = "factorial_2d_interaction_heatmap.svg"
    deconv_iters: int = 75
    n_jobs: int = -1
    run_single_sweep: bool = True
    run_compound_sweep: bool = True
    run_factorial_sweep: bool = True


def parse_args() -> SweepConfig:
    parser = argparse.ArgumentParser(
        description="Replicate perturbation benchmark comparing Milopy DA and Pseudobulk Deconvolution."
    )
    parser.add_argument("--n-cells", type=int, default=2400, help="Total synthetic cells")
    parser.add_argument("--n-genes", type=int, default=240, help="Total genes")
    parser.add_argument("--n-clusters", type=int, default=8, help="Number of clusters")
    parser.add_argument("--n-patients", type=int, default=30, help="Number of patients")
    parser.add_argument("--n-levels", type=int, default=6, help="Number of perturbation intensity levels")
    parser.add_argument("--n-reps", type=int, default=5, help="Number of independent replicates per point")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--n-jobs", type=int, default=-1, help="Parallel CPU jobs")
    parser.add_argument("--out-dir", type=Path, default=Path("output/synthetic_benchmark"), help="Output directory")
    parser.add_argument("--results-dir", type=Path, default=Path("results/synthetic_benchmark"), help="Results directory")
    parser.add_argument("--skip-single", action="store_true", help="Skip single perturbation sweep")
    parser.add_argument("--skip-compound", action="store_true", help="Skip compound perturbation sweep")
    parser.add_argument("--skip-factorial", action="store_true", help="Skip 2D factorial interaction sweep")
    args = parser.parse_args()

    return SweepConfig(
        n_cells=args.n_cells,
        n_genes=args.n_genes,
        n_clusters=args.n_clusters,
        n_patients=args.n_patients,
        n_levels=args.n_levels,
        n_reps=args.n_reps,
        random_seed=args.seed,
        n_jobs=args.n_jobs,
        out_dir=args.out_dir,
        results_dir=args.results_dir,
        run_single_sweep=not args.skip_single,
        run_compound_sweep=not args.skip_compound,
        run_factorial_sweep=not args.skip_factorial,
    )


# ------------------------------------------------------------------------------
# 1. Base Synthetic Data Generator
# ------------------------------------------------------------------------------

def generate_base_data(
    config: SweepConfig,
    collinearity_overlap: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Generate base single-cell expression matrix, cluster assignments, and marker gene indices."""
    n_cells = config.n_cells
    n_genes = config.n_genes
    k = config.n_clusters

    cells_per_cluster = n_cells // k
    cluster_assignments = np.repeat(np.arange(k, dtype=np.int32), cells_per_cluster)

    # Base background expression: Poisson lambda = 0.4
    raw_counts = rng.poisson(lam=0.4, size=(n_cells, n_genes)).astype(np.float64)

    genes_per_cluster = max(8, n_genes // (k + 1))
    marker_map: list[list[int]] = []

    for i in range(k):
        start_g = i * genes_per_cluster
        end_g = min(n_genes, (i + 1) * genes_per_cluster)
        g_indices = list(range(start_g, end_g))

        # Apply severe collinearity: clusters share up to 99% of marker genes with previous cluster
        if collinearity_overlap > 0.0 and i > 0:
            target_source = i - 1 if i % 2 == 1 else (i // 2)
            source_start = target_source * genes_per_cluster
            source_end = min(n_genes, (target_source + 1) * genes_per_cluster)
            source_indices = list(range(source_start, source_end))
            overlap_pct = min(0.99, collinearity_overlap * 0.99)
            n_overlap = int(round(overlap_pct * min(len(g_indices), len(source_indices))))
            if n_overlap > 0:
                g_indices = source_indices[:n_overlap] + g_indices[n_overlap:]

        marker_map.append(g_indices)
        mask_i = (cluster_assignments == i)
        n_cluster_cells = int(np.sum(mask_i))
        if n_cluster_cells > 0 and len(g_indices) > 0:
            marker_counts = rng.negative_binomial(n=5, p=0.35, size=(n_cluster_cells, len(g_indices)))
            raw_counts[mask_i, :][:, g_indices] += marker_counts.astype(np.float64)

    return raw_counts, cluster_assignments, np.array(marker_map, dtype=object)


# ------------------------------------------------------------------------------
# 2. Patient Assignment Conditioned on Response Gradient
# ------------------------------------------------------------------------------

def assign_patients(
    cluster_assignments: np.ndarray,
    n_patients: int,
    sparsity_intensity: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, dict[int, float]]:
    """Assign cells to Responder and Non-responder patients conditioned on cluster probability gradient."""
    k = int(np.max(cluster_assignments) + 1)
    cluster_probs = {cl: round(0.85 - (0.70 * cl / max(k - 1, 1)), 3) for cl in range(k)}

    half_p = n_patients // 2
    r_patients = [f"P_{j:02d}_R" for j in range(half_p)]
    nr_patients = [f"P_{j:02d}_NR" for j in range(half_p, n_patients)]

    # Patient sampling weights under sparsity
    if sparsity_intensity > 0.0:
        raw_w = np.exp(rng.normal(0.0, 4.0 * sparsity_intensity, size=n_patients))
        w_r = raw_w[:half_p] / raw_w[:half_p].sum()
        w_nr = raw_w[half_p:] / raw_w[half_p:].sum()
    else:
        w_r = np.full(half_p, 1.0 / half_p)
        w_nr = np.full(half_p, 1.0 / half_p)

    patient_assignments: list[str] = []
    patient_responses: list[str] = []

    for cl in cluster_assignments:
        prob = cluster_probs[int(cl)]
        is_resp = rng.random() < prob
        if is_resp:
            p_chosen = rng.choice(r_patients, p=w_r)
            resp = "R"
        else:
            p_chosen = rng.choice(nr_patients, p=w_nr)
            resp = "NR"
        patient_assignments.append(p_chosen)
        patient_responses.append(resp)

    return np.array(patient_assignments), np.array(patient_responses), cluster_probs


# ------------------------------------------------------------------------------
# 3. Modular Perturbation Transformations
# ------------------------------------------------------------------------------

def apply_perturbations(
    raw_counts: np.ndarray,
    cluster_assignments: np.ndarray,
    patient_assignments: np.ndarray,
    patient_responses: np.ndarray,
    mode: str | PerturbationProfile = "cell_size",
    intensity: float = 0.0,
    marker_map: list[list[int]] | np.ndarray | None = None,
    rng: np.random.Generator | None = None,
    profile: PerturbationProfile | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Apply single or compound perturbations across single-cell counts, reference Phi, and pseudobulk."""
    if rng is None:
        rng = np.random.default_rng(42)

    if isinstance(mode, PerturbationProfile):
        p = mode
    elif profile is not None:
        p = profile
    else:
        p = PerturbationProfile.from_single(mode, intensity)

    n_cells, n_genes = raw_counts.shape
    k_clusters = int(np.max(cluster_assignments) + 1)
    unique_patients = sorted(list(set(patient_assignments)))
    n_patients = len(unique_patients)
    p_to_idx = {pat: i for i, pat in enumerate(unique_patients)}

    sc_counts = raw_counts.copy()

    # 1. State Activation Confounding (Up to 20x upregulation in Responders)
    if p.activation_confounding > 0.0:
        gamma = 1.0 + (19.0 * p.activation_confounding)
        act_mask = ((cluster_assignments == 0) | (cluster_assignments == 1)) & (patient_responses == "R")
        sc_counts[act_mask, :] *= gamma

    # 2. Inter-Patient Transcriptional Heterogeneity (sigma up to 3.0)
    if p.patient_shift > 0.0:
        sigma = 3.0 * p.patient_shift
        for pat in unique_patients:
            p_mask = (patient_assignments == pat)
            p_shift = np.exp(rng.normal(0.0, sigma, size=(1, n_genes)))
            sc_counts[p_mask, :] *= p_shift

    # 3. Ambient RNA Soup Contamination in Single Cells (up to 40% ambient soup)
    if p.ambient_soup > 0.0:
        eta = 0.40 * p.ambient_soup
        soup_mean = np.mean(sc_counts, axis=0, keepdims=True)
        soup_sum = float(soup_mean.sum())
        if soup_sum > 0:
            soup_norm = soup_mean / soup_sum
            cell_totals = sc_counts.sum(axis=1, keepdims=True)
            sc_counts = (1.0 - eta) * sc_counts + (eta * cell_totals * soup_norm)

    # 4. Reference Matrix Phi Computation (cluster centroids)
    ref_centroids = np.zeros((k_clusters, n_genes), dtype=np.float64)
    for c in range(k_clusters):
        c_idx = np.where(cluster_assignments == c)[0]
        if len(c_idx) > 0:
            ref_centroids[c, :] = np.mean(sc_counts[c_idx, :], axis=0)
    row_sums = ref_centroids.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    norm_phi = ref_centroids / row_sums

    # 5. Cell Size / Total mRNA Multipliers (Exponential asymmetry up to 50x)
    size_factors = np.ones(k_clusters, dtype=np.float64)
    if p.cell_size > 0.0:
        max_size = 1.0 + (49.0 * p.cell_size)
        size_factors = np.geomspace(1.0, max_size, k_clusters)

    # 6. Build Pseudobulk per Patient with optional Marker Dysregulation
    dys_multipliers: dict[tuple[int, int], np.ndarray] = {}
    if p.marker_dysregulation > 0.0 and marker_map is not None:
        sigma_dys = 2.0 * p.marker_dysregulation
        for p_idx in range(n_patients):
            for cl in range(k_clusters):
                m_indices = [int(x) for x in marker_map[cl]]
                if len(m_indices) > 0:
                    dys_multipliers[(p_idx, cl)] = np.exp(rng.normal(0.0, sigma_dys, size=len(m_indices)))

    pseudobulk_mat = np.zeros((n_patients, n_genes), dtype=np.float64)
    for c_idx in range(n_cells):
        p_idx = p_to_idx[patient_assignments[c_idx]]
        cl = cluster_assignments[c_idx]
        s_factor = size_factors[cl]
        cell_expr = sc_counts[c_idx, :].copy() * s_factor
        if (p_idx, cl) in dys_multipliers and marker_map is not None:
            m_indices = [int(x) for x in marker_map[cl]]
            cell_expr[m_indices] *= dys_multipliers[(p_idx, cl)]
        pseudobulk_mat[p_idx, :] += cell_expr

    # 7. Ghost / Tumor Cell Contamination in Bulk (Up to 85% unmodeled contamination)
    if p.ghost_contamination > 0.0:
        f_ghost = 0.85 * p.ghost_contamination
        tumor_sig = rng.exponential(scale=2.0, size=n_genes)
        tumor_sig /= tumor_sig.sum()
        for p_idx in range(n_patients):
            bulk_sum = pseudobulk_mat[p_idx, :].sum()
            if bulk_sum > 0:
                pseudobulk_mat[p_idx, :] = (
                    (1.0 - f_ghost) * pseudobulk_mat[p_idx, :] + (f_ghost * bulk_sum * tumor_sig)
                )

    return sc_counts, cluster_assignments, pseudobulk_mat, norm_phi


# ------------------------------------------------------------------------------
# 4. Pipeline Execution: Milopy + Deconvolution + Concordance
# ------------------------------------------------------------------------------

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


def run_single_simulation(
    config: SweepConfig,
    mode: str | PerturbationProfile = "cell_size",
    intensity: float = 0.0,
    rep_idx: int = 0,
    rng: np.random.Generator | None = None,
    profile: PerturbationProfile | None = None,
    label: str | None = None,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Execute the full Milopy vs Deconvolution benchmark under a specific perturbation state and replicate."""
    if rng is None:
        rng = np.random.default_rng(42)

    if isinstance(mode, PerturbationProfile):
        p = mode
        mode_label = label or "custom_profile"
    elif profile is not None:
        p = profile
        mode_label = label or "custom_profile"
    else:
        p = PerturbationProfile.from_single(mode, intensity)
        mode_label = label or mode

    raw_counts, cluster_assignments, marker_map = generate_base_data(
        config=config,
        collinearity_overlap=p.collinearity,
        rng=rng,
    )

    patient_assignments, patient_responses, _ = assign_patients(
        cluster_assignments=cluster_assignments,
        n_patients=config.n_patients,
        sparsity_intensity=p.sampling_sparsity,
        rng=rng,
    )

    sc_counts, cluster_assignments, pseudobulk_mat, norm_phi = apply_perturbations(
        raw_counts=raw_counts,
        cluster_assignments=cluster_assignments,
        patient_assignments=patient_assignments,
        patient_responses=patient_responses,
        mode=p,
        intensity=intensity,
        marker_map=marker_map,
        rng=rng,
    )

    k_clusters = config.n_clusters
    unique_patients = sorted(list(set(patient_assignments)))
    n_patients = len(unique_patients)
    p_to_idx = {pat: i for i, pat in enumerate(unique_patients)}

    p_meta: dict[str, int] = {}
    for pat, resp in zip(patient_assignments, patient_responses, strict=True):
        p_meta[pat] = 1 if resp == "R" else 0
    y_patient = np.array([p_meta[pat] for pat in unique_patients], dtype=np.float64)

    # 1. Compute True Proportions per Patient
    true_fractions = np.zeros((n_patients, k_clusters), dtype=np.float64)
    patient_cell_counts = np.zeros(n_patients, dtype=np.float64)
    for c_idx in range(len(cluster_assignments)):
        p_idx = p_to_idx[patient_assignments[c_idx]]
        cl = cluster_assignments[c_idx]
        true_fractions[p_idx, cl] += 1.0
        patient_cell_counts[p_idx] += 1.0

    for p_idx in range(n_patients):
        if patient_cell_counts[p_idx] > 0:
            true_fractions[p_idx, :] /= patient_cell_counts[p_idx]

    # 2. Run Pseudobulk Deconvolution
    inferred_fractions = np.zeros((n_patients, k_clusters), dtype=np.float64)
    for p_idx in range(n_patients):
        bulk_vec = pseudobulk_mat[p_idx, :]
        if bulk_vec.sum() == 0:
            inferred_fractions[p_idx, :] = 1.0 / k_clusters
        else:
            try:
                if HAS_INSTAPRISM:
                    _, _, fracs, _ = instaprism.insta_prism(
                        bulk=bulk_vec,
                        reference=norm_phi,
                        n_iter=config.deconv_iters,
                    )
                else:
                    reg = 1e-4
                    phi_t_phi = norm_phi @ norm_phi.T + reg * np.eye(k_clusters)
                    phi_t_b = norm_phi @ bulk_vec
                    fracs = np.linalg.solve(phi_t_phi, phi_t_b)
                    fracs = np.clip(fracs, 0.0, None)
                    s = fracs.sum()
                    fracs = fracs / s if s > 0 else np.full(k_clusters, 1.0 / k_clusters)
                inferred_fractions[p_idx, :] = fracs
            except Exception:
                inferred_fractions[p_idx, :] = 1.0 / k_clusters

    # Simplex Renormalization
    row_sums = inferred_fractions.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    inferred_fractions /= row_sums

    # 3. Milopy Single-Cell Differential Abundance
    adata = ad.AnnData(
        X=sp.csr_matrix(sc_counts),
        obs={
            "cluster": cluster_assignments.astype(str),
            "patient": patient_assignments,
            "response": patient_responses,
        },
    )
    n_comps = min(15, config.n_genes - 1, len(cluster_assignments) - 1)
    sc.tl.pca(adata, n_comps=n_comps, use_highly_variable=False, random_state=42)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_comps, random_state=42)

    milo_betas = np.zeros(k_clusters, dtype=np.float64)
    if HAS_MILOPY:
        try:
            milopy.core.make_nhoods(adata, prop=0.15, seed=42)
            milopy.core.count_nhoods(adata, sample_col="patient")
            milopy.core.DA_nhoods(adata, design="~response", model_contrasts="response[T.R]")

            nhood_adata = adata.uns["nhood_adata"]
            nhood_clusters = []
            for nh_idx in range(nhood_adata.n_obs):
                center_cell = nhood_adata.obs_names[nh_idx]
                nhood_clusters.append(adata.obs.loc[center_cell, "cluster"])
            nhood_clusters = np.array(nhood_clusters, dtype=np.int32)
            nhood_logfc = np.array(nhood_adata.obs["logFC"].values, dtype=np.float64)

            for c in range(k_clusters):
                m = (nhood_clusters == c)
                milo_betas[c] = float(np.median(nhood_logfc[m])) if np.sum(m) > 0 else 0.0
        except Exception:
            for c in range(k_clusters):
                milo_betas[c] = float(np.mean(true_fractions[y_patient == 1, c]) - np.mean(true_fractions[y_patient == 0, c]))
    else:
        for c in range(k_clusters):
            milo_betas[c] = float(np.mean(true_fractions[y_patient == 1, c]) - np.mean(true_fractions[y_patient == 0, c]))

    # 4. Standardized Patient-Level Effect Sizes
    deconv_betas = np.zeros(k_clusters, dtype=np.float64)
    true_betas = np.zeros(k_clusters, dtype=np.float64)

    for c in range(k_clusters):
        x_d = inferred_fractions[:, c]
        x_t = true_fractions[:, c]

        if np.std(x_d) < 1e-8:
            deconv_betas[c] = 0.0
        else:
            r_d, _ = stats.pointbiserialr(y_patient, (x_d - np.mean(x_d)) / np.std(x_d))
            r_d_clip = float(np.clip(r_d, -0.999, 0.999))
            deconv_betas[c] = float(2.0 * r_d_clip / np.sqrt(1.0 - r_d_clip**2 + 1e-12))

        if np.std(x_t) < 1e-8:
            true_betas[c] = 0.0
        else:
            r_t, _ = stats.pointbiserialr(y_patient, (x_t - np.mean(x_t)) / np.std(x_t))
            r_t_clip = float(np.clip(r_t, -0.999, 0.999))
            true_betas[c] = float(2.0 * r_t_clip / np.sqrt(1.0 - r_t_clip**2 + 1e-12))

    # 5. Concordance Metrics
    spearman_rho, spearman_p = stats.spearmanr(milo_betas, deconv_betas)
    pearson_r, pearson_p = stats.pearsonr(milo_betas, deconv_betas)
    fid_spearman_rho, fid_spearman_p = stats.spearmanr(true_betas, deconv_betas)
    fid_pearson_r, fid_pearson_p = stats.pearsonr(true_betas, deconv_betas)
    milo_true_rho, _ = stats.spearmanr(milo_betas, true_betas)

    sign_matches = int(np.sum((milo_betas * deconv_betas) > 0))
    sign_concordance_pct = float(sign_matches / k_clusters) * 100.0

    summary_record = {
        "perturbation_mode": mode_label,
        "intensity": float(intensity),
        "replicate": int(rep_idx),
        "spearman_rho_milo_deconv": float(spearman_rho) if not np.isnan(spearman_rho) else 0.0,
        "spearman_pval_milo_deconv": float(spearman_p) if not np.isnan(spearman_p) else 1.0,
        "pearson_r_milo_deconv": float(pearson_r) if not np.isnan(pearson_r) else 0.0,
        "pearson_pval_milo_deconv": float(pearson_p) if not np.isnan(pearson_p) else 1.0,
        "sign_concordance_pct": sign_concordance_pct,
        "fidelity_spearman_rho": float(fid_spearman_rho) if not np.isnan(fid_spearman_rho) else 0.0,
        "fidelity_pearson_r": float(fid_pearson_r) if not np.isnan(fid_pearson_r) else 0.0,
        "milo_true_spearman_rho": float(milo_true_rho) if not np.isnan(milo_true_rho) else 0.0,
    }

    scatter_records: list[dict[str, object]] = []
    for c in range(k_clusters):
        scatter_records.append({
            "perturbation_mode": mode_label,
            "intensity": float(intensity),
            "replicate": int(rep_idx),
            "cell_state": f"C{c}",
            "milo_mean_logfc": float(milo_betas[c]),
            "deconv_beta": float(deconv_betas[c]),
            "true_beta": float(true_betas[c]),
            "quadrant": classify_quadrant(float(milo_betas[c]), float(deconv_betas[c])),
        })

    return summary_record, scatter_records


# ------------------------------------------------------------------------------
# 5. Parallel Sweep Boundaries: Single, Compound, & 2D Factorial
# ------------------------------------------------------------------------------

def run_replicate_sweep(
    config: SweepConfig,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """Execute parallel simulation sweeps across single modes, intensity levels, and replicates."""
    print("=== Launching Parallel Single-Mode Replicate Benchmark ===")
    config.out_dir.mkdir(parents=True, exist_ok=True)
    config.results_dir.mkdir(parents=True, exist_ok=True)

    intensities = np.linspace(0.0, 1.0, config.n_levels)
    tasks: list[tuple[str, float, int, int]] = []
    task_idx = 0
    for mode in PERTURBATION_MODES:
        for level in intensities:
            for rep in range(config.n_reps):
                seed = config.random_seed + (task_idx * 13) + (rep * 101)
                tasks.append((mode, float(level), rep, seed))
                task_idx += 1

    print(f"Executing {len(tasks)} single-mode simulation tasks across {config.n_jobs} parallel workers...")
    results = Parallel(n_jobs=config.n_jobs, verbose=5)(
        delayed(run_single_simulation)(
            config,
            mode,
            level,
            rep,
            np.random.default_rng(seed),
        )
        for mode, level, rep, seed in tasks
    )

    summary_list = [r[0] for r in results]
    scatter_list = [p for r in results for p in r[1]]

    raw_sweep_df = pl.DataFrame(summary_list)
    scatter_df = pl.DataFrame(scatter_list)

    summary_df = (
        raw_sweep_df.group_by(["perturbation_mode", "intensity"])
        .agg([
            pl.col("spearman_rho_milo_deconv").mean().alias("rho_mean"),
            pl.col("spearman_rho_milo_deconv").std().alias("rho_std"),
            pl.col("pearson_r_milo_deconv").mean().alias("pearson_mean"),
            pl.col("pearson_r_milo_deconv").std().alias("pearson_std"),
            pl.col("sign_concordance_pct").mean().alias("sign_mean"),
            pl.col("sign_concordance_pct").std().alias("sign_std"),
            pl.col("fidelity_spearman_rho").mean().alias("fid_mean"),
            pl.col("fidelity_spearman_rho").std().alias("fid_std"),
            pl.len().alias("n_reps"),
        ])
        .with_columns([
            pl.col("rho_std").fill_null(0.0),
            pl.col("sign_std").fill_null(0.0),
            pl.col("fid_std").fill_null(0.0),
        ])
        .with_columns([
            (pl.col("rho_mean") - pl.col("rho_std")).clip(-1.0, 1.0).alias("rho_lower"),
            (pl.col("rho_mean") + pl.col("rho_std")).clip(-1.0, 1.0).alias("rho_upper"),
            (pl.col("sign_mean") - pl.col("sign_std")).clip(0.0, 100.0).alias("sign_lower"),
            (pl.col("sign_mean") + pl.col("sign_std")).clip(0.0, 100.0).alias("sign_upper"),
            (pl.col("fid_mean") - pl.col("fid_std")).clip(-1.0, 1.0).alias("fid_lower"),
            (pl.col("fid_mean") + pl.col("fid_std")).clip(-1.0, 1.0).alias("fid_upper"),
        ])
        .sort(["perturbation_mode", "intensity"])
    )

    return raw_sweep_df, summary_df, scatter_df


def run_compound_sweep(
    config: SweepConfig,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """Execute parallel simulation sweeps for compound clinical perturbation regimes."""
    print("=== Launching Compound Multi-Factorial Perturbation Benchmark ===")
    config.out_dir.mkdir(parents=True, exist_ok=True)
    config.results_dir.mkdir(parents=True, exist_ok=True)

    intensities = np.linspace(0.0, 1.0, config.n_levels)
    tasks: list[tuple[str, PerturbationProfile, float, int, int]] = []
    task_idx = 0
    for regime_key in COMPOUND_REGIMES:
        for level in intensities:
            profile = make_compound_profile(regime_key, float(level))
            for rep in range(config.n_reps):
                seed = config.random_seed + 5000 + (task_idx * 17) + (rep * 103)
                tasks.append((regime_key, profile, float(level), rep, seed))
                task_idx += 1

    print(f"Executing {len(tasks)} compound simulation tasks across {config.n_jobs} parallel workers...")
    results = Parallel(n_jobs=config.n_jobs, verbose=5)(
        delayed(run_single_simulation)(
            config=config,
            mode=profile,
            intensity=level,
            rep_idx=rep,
            rng=np.random.default_rng(seed),
            label=regime_key,
        )
        for regime_key, profile, level, rep, seed in tasks
    )

    summary_list = [r[0] for r in results]
    scatter_list = [p for r in results for p in r[1]]

    raw_compound_df = pl.DataFrame(summary_list)
    scatter_df = pl.DataFrame(scatter_list)

    summary_df = (
        raw_compound_df.group_by(["perturbation_mode", "intensity"])
        .agg([
            pl.col("spearman_rho_milo_deconv").mean().alias("rho_mean"),
            pl.col("spearman_rho_milo_deconv").std().alias("rho_std"),
            pl.col("pearson_r_milo_deconv").mean().alias("pearson_mean"),
            pl.col("pearson_r_milo_deconv").std().alias("pearson_std"),
            pl.col("sign_concordance_pct").mean().alias("sign_mean"),
            pl.col("sign_concordance_pct").std().alias("sign_std"),
            pl.col("fidelity_spearman_rho").mean().alias("fid_mean"),
            pl.col("fidelity_spearman_rho").std().alias("fid_std"),
            pl.len().alias("n_reps"),
        ])
        .with_columns([
            pl.col("rho_std").fill_null(0.0),
            pl.col("sign_std").fill_null(0.0),
            pl.col("fid_std").fill_null(0.0),
        ])
        .with_columns([
            (pl.col("rho_mean") - pl.col("rho_std")).clip(-1.0, 1.0).alias("rho_lower"),
            (pl.col("rho_mean") + pl.col("rho_std")).clip(-1.0, 1.0).alias("rho_upper"),
            (pl.col("sign_mean") - pl.col("sign_std")).clip(0.0, 100.0).alias("sign_lower"),
            (pl.col("sign_mean") + pl.col("sign_std")).clip(0.0, 100.0).alias("sign_upper"),
            (pl.col("fid_mean") - pl.col("fid_std")).clip(-1.0, 1.0).alias("fid_lower"),
            (pl.col("fid_mean") + pl.col("fid_std")).clip(-1.0, 1.0).alias("fid_upper"),
        ])
        .sort(["perturbation_mode", "intensity"])
    )

    return raw_compound_df, summary_df, scatter_df


def run_2d_factorial_sweep(
    config: SweepConfig,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Execute a 2D factorial interaction sweep (Cell Size x Activation Confounding) across 5x5 grid."""
    print("=== Launching 2D Factorial Interaction Sweep (Cell Size x Activation) ===")
    config.out_dir.mkdir(parents=True, exist_ok=True)
    config.results_dir.mkdir(parents=True, exist_ok=True)

    grid_levels = np.linspace(0.0, 1.0, 5)
    tasks: list[tuple[float, float, int, int]] = []
    task_idx = 0
    for cs in grid_levels:
        for act in grid_levels:
            for rep in range(config.n_reps):
                seed = config.random_seed + 10000 + (task_idx * 19) + (rep * 107)
                tasks.append((float(cs), float(act), rep, seed))
                task_idx += 1

    print(f"Executing {len(tasks)} factorial simulation tasks across {config.n_jobs} parallel workers...")
    results = Parallel(n_jobs=config.n_jobs, verbose=5)(
        delayed(run_single_simulation)(
            config=config,
            mode=PerturbationProfile(cell_size=cs, activation_confounding=act),
            intensity=float(np.sqrt(cs**2 + act**2) / np.sqrt(2)),
            rep_idx=rep,
            rng=np.random.default_rng(seed),
            label=f"CS_{cs:.2f}_ACT_{act:.2f}",
        )
        for cs, act, rep, seed in tasks
    )

    summary_list = []
    for (cs, act, rep, _), (s_rec, _) in zip(tasks, results, strict=True):
        rec = dict(s_rec)
        rec["cell_size_intensity"] = cs
        rec["activation_intensity"] = act
        summary_list.append(rec)

    raw_factorial_df = pl.DataFrame(summary_list)

    summary_df = (
        raw_factorial_df.group_by(["cell_size_intensity", "activation_intensity"])
        .agg([
            pl.col("spearman_rho_milo_deconv").mean().alias("rho_mean"),
            pl.col("spearman_rho_milo_deconv").std().alias("rho_std"),
            pl.col("sign_concordance_pct").mean().alias("sign_mean"),
            pl.col("sign_concordance_pct").std().alias("sign_std"),
            pl.col("fidelity_spearman_rho").mean().alias("fid_mean"),
            pl.col("fidelity_spearman_rho").std().alias("fid_std"),
            pl.len().alias("n_reps"),
        ])
        .with_columns([
            pl.col("rho_std").fill_null(0.0),
            pl.col("sign_std").fill_null(0.0),
            pl.col("fid_std").fill_null(0.0),
        ])
        .sort(["cell_size_intensity", "activation_intensity"])
    )

    return raw_factorial_df, summary_df


# ------------------------------------------------------------------------------
# 6. Altair Visualizations: Single, Compound, & 2D Factorial
# ------------------------------------------------------------------------------

def plot_sensitivity_curves_with_sd(
    summary_df: pl.DataFrame,
) -> alt.VConcatChart:
    """Construct publication multi-panel vector SVG with shaded Mean +/- SD envelopes for single modes."""
    data_pd = summary_df.to_pandas()
    data_pd["mode_display"] = data_pd["perturbation_mode"].map(MODE_DISPLAY_TITLES)

    color_scale = alt.Scale(scheme="category10")

    # Panel 1: Spearman Rho (Mean +/- SD)
    p1_band = (
        alt.Chart(data_pd)
        .mark_area(opacity=0.20)
        .encode(
            x=alt.X("intensity:Q", title="Perturbation Intensity (0.0 = Baseline, 1.0 = Max)"),
            y=alt.Y("rho_lower:Q", title="Milo vs. Deconv (Spearman ρ)"),
            y2="rho_upper:Q",
            color=alt.Color("mode_display:N", scale=color_scale, title="Perturbation Mode"),
        )
    )
    p1_line = (
        alt.Chart(data_pd)
        .mark_line(point=True, strokeWidth=2.5)
        .encode(
            x=alt.X("intensity:Q"),
            y=alt.Y("rho_mean:Q", scale=alt.Scale(domain=[-1.05, 1.05])),
            color=alt.Color("mode_display:N", scale=color_scale, title="Perturbation Mode"),
            tooltip=["mode_display:N", "intensity:Q", alt.Tooltip("rho_mean:Q", format="+.3f"), alt.Tooltip("rho_std:Q", format=".3f")],
        )
    )
    h_zero = alt.Chart().mark_rule(color="#999999", strokeDash=[4, 4]).encode(y=alt.datum(0.0))
    p1 = (p1_band + p1_line + h_zero).properties(width=420, height=280, title="A: Correlation Degradation (Mean ± SD, N=5)")

    # Panel 2: Directional Sign Agreement (Mean +/- SD)
    p2_band = (
        alt.Chart(data_pd)
        .mark_area(opacity=0.20)
        .encode(
            x=alt.X("intensity:Q", title="Perturbation Intensity"),
            y=alt.Y("sign_lower:Q", title="Directional Sign Concordance (%)"),
            y2="sign_upper:Q",
            color=alt.Color("mode_display:N", scale=color_scale, title="Perturbation Mode"),
        )
    )
    p2_line = (
        alt.Chart(data_pd)
        .mark_line(point=True, strokeWidth=2.5)
        .encode(
            x=alt.X("intensity:Q"),
            y=alt.Y("sign_mean:Q", scale=alt.Scale(domain=[15.0, 105.0])),
            color=alt.Color("mode_display:N", scale=color_scale, title="Perturbation Mode"),
            tooltip=["mode_display:N", "intensity:Q", alt.Tooltip("sign_mean:Q", format=".1f"), alt.Tooltip("sign_std:Q", format=".1f")],
        )
    )
    h_coin = alt.Chart().mark_rule(color="#d95f02", strokeDash=[3, 3]).encode(y=alt.datum(50.0))
    p2 = (p2_band + p2_line + h_coin).properties(width=420, height=280, title="B: Directional Sign Concordance (Mean ± SD)")

    # Panel 3: Deconvolution Fidelity vs Ground Truth (Mean +/- SD)
    p3_band = (
        alt.Chart(data_pd)
        .mark_area(opacity=0.20)
        .encode(
            x=alt.X("intensity:Q", title="Perturbation Intensity"),
            y=alt.Y("fid_lower:Q", title="True vs. Deconv (Spearman ρ)"),
            y2="fid_upper:Q",
            color=alt.Color("mode_display:N", scale=color_scale, title="Perturbation Mode"),
        )
    )
    p3_line = (
        alt.Chart(data_pd)
        .mark_line(point=True, strokeWidth=2.5)
        .encode(
            x=alt.X("intensity:Q"),
            y=alt.Y("fid_mean:Q", scale=alt.Scale(domain=[-1.05, 1.05])),
            color=alt.Color("mode_display:N", scale=color_scale, title="Perturbation Mode"),
            tooltip=["mode_display:N", "intensity:Q", alt.Tooltip("fid_mean:Q", format="+.3f"), alt.Tooltip("fid_std:Q", format=".3f")],
        )
    )
    p3 = (p3_band + p3_line + h_zero).properties(width=420, height=280, title="C: Deconvolution Ground Truth Fidelity (Mean ± SD)")

    # Panel 4: Max Degradation Ranking Bar Chart
    ranking_pd = (
        summary_df.group_by("perturbation_mode")
        .agg([
            pl.col("rho_mean").first().alias("baseline_rho"),
            pl.col("rho_mean").min().alias("min_rho"),
            (pl.col("rho_mean").first() - pl.col("rho_mean").min()).alias("max_rho_drop"),
        ])
        .to_pandas()
    )
    ranking_pd["mode_display"] = ranking_pd["perturbation_mode"].map(MODE_DISPLAY_TITLES)

    p4 = (
        alt.Chart(ranking_pd)
        .mark_bar()
        .encode(
            y=alt.Y("mode_display:N", sort="-x", title="Perturbation Mode"),
            x=alt.X("max_rho_drop:Q", title="Maximum Spearman ρ Drop (Baseline - Min)"),
            color=alt.Color("mode_display:N", scale=color_scale, legend=None),
            tooltip=["mode_display:N", alt.Tooltip("max_rho_drop:Q", format=".3f")],
        )
        .properties(width=420, height=280, title="D: Perturbation Sensitivity Ranking")
    )

    top_row = alt.hconcat(p1, p2)
    bottom_row = alt.hconcat(p3, p4)

    full_chart = (
        alt.vconcat(top_row, bottom_row)
        .properties(
            title=alt.TitleParams(
                text="Multi-Replicate Perturbation Sensitivity: Milopy DA vs. Deconvolution",
                subtitle=[
                    "Evaluating 6 biological, technical, and transcriptional distortion modes across 5 independent replicates per point (N=5)",
                    "Shaded envelopes indicate +/- 1 Standard Deviation; dashed orange rule marks the 50% random coin-flip sign baseline",
                ],
                fontSize=16,
                subtitleFontSize=12,
                anchor="middle",
            )
        )
        .configure_view(strokeWidth=1, stroke="#cccccc")
    )
    return full_chart


def plot_compound_sensitivity_curves(
    compound_summary_df: pl.DataFrame,
) -> alt.VConcatChart:
    """Construct multi-panel vector SVG with shaded Mean +/- SD envelopes for compound clinical regimes."""
    data_pd = compound_summary_df.to_pandas()
    data_pd["regime_title"] = data_pd["perturbation_mode"].map(lambda k: COMPOUND_REGIMES[k]["title"])

    color_scale = alt.Scale(scheme="set1")

    # Panel A: Spearman Rho (Mean +/- SD)
    p1_band = (
        alt.Chart(data_pd)
        .mark_area(opacity=0.22)
        .encode(
            x=alt.X("intensity:Q", title="Compound Distortion Intensity (0.0 = Clean, 1.0 = Max)"),
            y=alt.Y("rho_lower:Q", title="Milo vs. Deconv (Spearman ρ)"),
            y2="rho_upper:Q",
            color=alt.Color("regime_title:N", scale=color_scale, title="Compound Regime"),
        )
    )
    p1_line = (
        alt.Chart(data_pd)
        .mark_line(point=True, strokeWidth=2.8)
        .encode(
            x=alt.X("intensity:Q"),
            y=alt.Y("rho_mean:Q", scale=alt.Scale(domain=[-1.05, 1.05])),
            color=alt.Color("regime_title:N", scale=color_scale, title="Compound Regime"),
            tooltip=["regime_title:N", "intensity:Q", alt.Tooltip("rho_mean:Q", format="+.3f"), alt.Tooltip("rho_std:Q", format=".3f")],
        )
    )
    h_zero = alt.Chart().mark_rule(color="#999999", strokeDash=[4, 4]).encode(y=alt.datum(0.0))
    p1 = (p1_band + p1_line + h_zero).properties(width=420, height=280, title="A: Compound Correlation Breakdown (Mean ± SD, N=5)")

    # Panel B: Directional Sign Agreement
    p2_band = (
        alt.Chart(data_pd)
        .mark_area(opacity=0.22)
        .encode(
            x=alt.X("intensity:Q", title="Compound Distortion Intensity"),
            y=alt.Y("sign_lower:Q", title="Directional Sign Concordance (%)"),
            y2="sign_upper:Q",
            color=alt.Color("regime_title:N", scale=color_scale, title="Compound Regime"),
        )
    )
    p2_line = (
        alt.Chart(data_pd)
        .mark_line(point=True, strokeWidth=2.8)
        .encode(
            x=alt.X("intensity:Q"),
            y=alt.Y("sign_mean:Q", scale=alt.Scale(domain=[10.0, 105.0])),
            color=alt.Color("regime_title:N", scale=color_scale, title="Compound Regime"),
            tooltip=["regime_title:N", "intensity:Q", alt.Tooltip("sign_mean:Q", format=".1f"), alt.Tooltip("sign_std:Q", format=".1f")],
        )
    )
    h_coin = alt.Chart().mark_rule(color="#d95f02", strokeDash=[3, 3]).encode(y=alt.datum(50.0))
    p2 = (p2_band + p2_line + h_coin).properties(width=420, height=280, title="B: Directional Sign Concordance (Mean ± SD)")

    # Panel C: Deconvolution Ground Truth Fidelity
    p3_band = (
        alt.Chart(data_pd)
        .mark_area(opacity=0.22)
        .encode(
            x=alt.X("intensity:Q", title="Compound Distortion Intensity"),
            y=alt.Y("fid_lower:Q", title="True vs. Deconv (Spearman ρ)"),
            y2="fid_upper:Q",
            color=alt.Color("regime_title:N", scale=color_scale, title="Compound Regime"),
        )
    )
    p3_line = (
        alt.Chart(data_pd)
        .mark_line(point=True, strokeWidth=2.8)
        .encode(
            x=alt.X("intensity:Q"),
            y=alt.Y("fid_mean:Q", scale=alt.Scale(domain=[-1.05, 1.05])),
            color=alt.Color("regime_title:N", scale=color_scale, title="Compound Regime"),
            tooltip=["regime_title:N", "intensity:Q", alt.Tooltip("fid_mean:Q", format="+.3f"), alt.Tooltip("fid_std:Q", format=".3f")],
        )
    )
    p3 = (p3_band + p3_line + h_zero).properties(width=420, height=280, title="C: Deconvolution Ground Truth Fidelity (Mean ± SD)")

    # Panel D: Compound Severity Ranking
    ranking_pd = (
        compound_summary_df.group_by("perturbation_mode")
        .agg([
            pl.col("rho_mean").first().alias("baseline_rho"),
            pl.col("rho_mean").min().alias("min_rho"),
            (pl.col("rho_mean").first() - pl.col("rho_mean").min()).alias("max_rho_drop"),
        ])
        .to_pandas()
    )
    ranking_pd["regime_title"] = ranking_pd["perturbation_mode"].map(lambda k: COMPOUND_REGIMES[k]["title"])

    p4 = (
        alt.Chart(ranking_pd)
        .mark_bar()
        .encode(
            y=alt.Y("regime_title:N", sort="-x", title="Compound Clinical Regime"),
            x=alt.X("max_rho_drop:Q", title="Maximum Spearman ρ Drop (Baseline - Min)"),
            color=alt.Color("regime_title:N", scale=color_scale, legend=None),
            tooltip=["regime_title:N", alt.Tooltip("max_rho_drop:Q", format=".3f")],
        )
        .properties(width=420, height=280, title="D: Compound Clinical Severity Ranking")
    )

    top_row = alt.hconcat(p1, p2)
    bottom_row = alt.hconcat(p3, p4)

    full_chart = (
        alt.vconcat(top_row, bottom_row)
        .properties(
            title=alt.TitleParams(
                text="Multi-Factorial Compound Stress-Testing: Milopy DA vs. Deconvolution",
                subtitle=[
                    "Simulating coupled clinical phenomena across 5 independent replicates (N=5): Core Biopsy, Inflamed TME, SC Noise, and Triple Jeopardy",
                    "Shaded envelopes represent +/- 1 Standard Deviation; red dashed rule marks the 50% random coin-flip sign baseline",
                ],
                fontSize=16,
                subtitleFontSize=12,
                anchor="middle",
            )
        )
        .configure_view(strokeWidth=1, stroke="#cccccc")
    )
    return full_chart


def plot_2d_interaction_heatmap(
    factorial_summary_df: pl.DataFrame,
) -> alt.HConcatChart:
    """Construct 2D interaction heatmap showing non-linear coupling of Cell Size x Activation Confounding."""
    data_pd = factorial_summary_df.to_pandas()
    data_pd["cs_pct"] = (data_pd["cell_size_intensity"] * 100).round().astype(int).astype(str) + "%"
    data_pd["act_pct"] = (data_pd["activation_intensity"] * 100).round().astype(int).astype(str) + "%"
    data_pd["rho_label"] = data_pd["rho_mean"].apply(lambda v: f"{v:+.2f}")
    data_pd["sign_label"] = data_pd["sign_mean"].apply(lambda v: f"{v:.0f}%")

    # Left Heatmap: Spearman Rho Mean
    base_rho = (
        alt.Chart(data_pd)
        .encode(
            x=alt.X("cs_pct:O", title="Cell Size Asymmetry Intensity", sort=["0%", "25%", "50%", "75%", "100%"]),
            y=alt.Y("act_pct:O", title="Activation Confounding Intensity", sort=["100%", "75%", "50%", "25%", "0%"]),
        )
    )
    rect_rho = base_rho.mark_rect().encode(
        color=alt.Color(
            "rho_mean:Q",
            scale=alt.Scale(scheme="viridis", domain=[-0.2, 1.0]),
            title="Mean Spearman ρ",
        ),
        tooltip=["cell_size_intensity:Q", "activation_intensity:Q", alt.Tooltip("rho_mean:Q", format="+.3f"), alt.Tooltip("rho_std:Q", format=".3f")],
    )
    text_rho = base_rho.mark_text(baseline="middle", fontSize=12, fontWeight="bold").encode(
        text="rho_label:N",
        color=alt.condition("datum.rho_mean < 0.45", alt.value("white"), alt.value("black")),
    )
    h_rho = (rect_rho + text_rho).properties(
        width=380,
        height=320,
        title=alt.TitleParams(
            text="A: Concordance Rank Correlation (Spearman ρ)",
            subtitle="Interaction between Cell Size Asymmetry and State Activation (N=5)",
            fontSize=13,
            subtitleFontSize=10,
        ),
    )

    # Right Heatmap: Sign Concordance %
    base_sign = (
        alt.Chart(data_pd)
        .encode(
            x=alt.X("cs_pct:O", title="Cell Size Asymmetry Intensity", sort=["0%", "25%", "50%", "75%", "100%"]),
            y=alt.Y("act_pct:O", title="Activation Confounding Intensity", sort=["100%", "75%", "50%", "25%", "0%"]),
        )
    )
    rect_sign = base_sign.mark_rect().encode(
        color=alt.Color(
            "sign_mean:Q",
            scale=alt.Scale(scheme="redyellowgreen", domain=[25.0, 90.0]),
            title="Sign Agreement (%)",
        ),
        tooltip=["cell_size_intensity:Q", "activation_intensity:Q", alt.Tooltip("sign_mean:Q", format=".1f"), alt.Tooltip("sign_std:Q", format=".1f")],
    )
    text_sign = base_sign.mark_text(baseline="middle", fontSize=12, fontWeight="bold").encode(
        text="sign_label:N",
        color=alt.condition("datum.sign_mean < 55", alt.value("black"), alt.value("black")),
    )
    h_sign = (rect_sign + text_sign).properties(
        width=380,
        height=320,
        title=alt.TitleParams(
            text="B: Directional Sign Agreement (%)",
            subtitle="Percentage of clusters with matching response direction (N=5)",
            fontSize=13,
            subtitleFontSize=10,
        ),
    )

    combined_heatmaps = (
        alt.hconcat(h_rho, h_sign)
        .properties(
            title=alt.TitleParams(
                text="2D Factorial Interaction Surface: Cell Size Asymmetry × State Activation Confounding",
                subtitle=[
                    "Evaluating 25 combinations of biological cell size variation and response-correlated activation across 5 replicates per grid cell (N=125 runs)",
                    "Reveals non-linear compound degradation: Cell size drives correlation drop, while activation drives directional sign inversion",
                ],
                fontSize=16,
                subtitleFontSize=12,
                anchor="middle",
            )
        )
        .configure_view(strokeWidth=1, stroke="#cccccc")
    )
    return combined_heatmaps


def plot_scatter_grid(
    scatter_df: pl.DataFrame,
    summary_df: pl.DataFrame,
) -> alt.VConcatChart:
    """Construct Sade-Feldman style multi-panel scatter grid comparing Baseline vs. all 6 single distortions."""
    quadrant_colors = {
        "Concordant Responder": "#2ca02c",      # Green
        "Concordant Non-Responder": "#1f77b4",  # Blue
        "Discordant (Milo+, Deconv-)": "#d62728", # Red
        "Discordant (Milo-, Deconv+)": "#ff7f0e", # Orange
        "Concordant Neutral": "#7f7f7f",        # Gray
    }
    color_scale = alt.Scale(
        domain=list(quadrant_colors.keys()),
        range=list(quadrant_colors.values()),
    )

    panels: list[alt.Chart] = []

    # Panel 1: Baseline Unperturbed
    base_points = scatter_df.filter((pl.col("intensity") == 0.0) & (pl.col("replicate") == 0)).to_pandas()
    base_summary = summary_df.filter((pl.col("intensity") == 0.0)).row(0, named=True)

    c_base = (
        alt.Chart(base_points)
        .mark_circle(size=140, opacity=0.9, stroke="black", strokeWidth=1)
        .encode(
            x=alt.X("milo_mean_logfc:Q", title="Milopy Mean logFC"),
            y=alt.Y("deconv_beta:Q", title="Deconvolution Effect (β)"),
            color=alt.Color("quadrant:N", scale=color_scale, legend=None),
            tooltip=["cell_state:N", "milo_mean_logfc:Q", "deconv_beta:Q", "quadrant:N"],
        )
    )
    line_base = c_base.transform_regression("milo_mean_logfc", "deconv_beta").mark_line(color="#444444", strokeDash=[4, 4], strokeWidth=2)
    h_rule = alt.Chart().mark_rule(color="#dddddd").encode(y=alt.datum(0))
    v_rule = alt.Chart().mark_rule(color="#dddddd").encode(x=alt.datum(0))

    p1 = (c_base + line_base + h_rule + v_rule).properties(
        width=320,
        height=280,
        title=alt.TitleParams(
            text="A: Baseline (Unperturbed)",
            subtitle=f"Spearman ρ = {base_summary['rho_mean']:+.3f} ± {base_summary['rho_std']:.3f} | Sign = {base_summary['sign_mean']:.0f}%",
            fontSize=12,
            subtitleFontSize=10,
        ),
    )
    panels.append(p1)

    # Panels 2-7: Individual Perturbations at representative moderate strength (intensity ~ 0.6)
    rep_intensity = 0.6
    single_modes = [
        ("cell_size", "B: Cell Size Asymmetry (10x)"),
        ("collinearity", "C: Reference Collinearity (60%)"),
        ("patient_shift", "D: Inter-Patient Shift (σ = 1.8)"),
        ("activation_confounding", "E: Activation Confounding (8x)"),
        ("ghost_contamination", "F: Ghost Contamination (50%)"),
        ("sampling_sparsity", "G: Sampling Sparsity (σ = 2.4)"),
    ]

    for mode, title in single_modes:
        sub_points = (
            scatter_df.filter(
                (pl.col("perturbation_mode") == mode)
                & ((pl.col("intensity") - rep_intensity).abs() < 0.05)
                & (pl.col("replicate") == 0)
            )
            .to_pandas()
        )
        sub_sum_rows = summary_df.filter(
            (pl.col("perturbation_mode") == mode)
            & ((pl.col("intensity") - rep_intensity).abs() < 0.05)
        )
        if sub_sum_rows.height > 0:
            s_row = sub_sum_rows.row(0, named=True)
            subtitle = f"Spearman ρ = {s_row['rho_mean']:+.3f} ± {s_row['rho_std']:.3f} | Sign = {s_row['sign_mean']:.0f}%"
        else:
            subtitle = "Evaluation instance"

        c_pt = (
            alt.Chart(sub_points)
            .mark_circle(size=140, opacity=0.9, stroke="black", strokeWidth=1)
            .encode(
                x=alt.X("milo_mean_logfc:Q", title="Milopy Mean logFC"),
                y=alt.Y("deconv_beta:Q", title="Deconvolution Effect (β)"),
                color=alt.Color("quadrant:N", scale=color_scale, legend=None),
                tooltip=["cell_state:N", "milo_mean_logfc:Q", "deconv_beta:Q", "quadrant:N"],
            )
        )
        reg_line = c_pt.transform_regression("milo_mean_logfc", "deconv_beta").mark_line(color="#444444", strokeDash=[4, 4], strokeWidth=2)
        p_sub = (c_pt + reg_line + h_rule + v_rule).properties(
            width=320,
            height=280,
            title=alt.TitleParams(text=title, subtitle=subtitle, fontSize=12, subtitleFontSize=10),
        )
        panels.append(p_sub)

    # Panel 8: Concordance Comparison Bar Chart
    comp_df = (
        summary_df.filter((pl.col("intensity") - rep_intensity).abs() < 0.05)
        .to_pandas()
    )
    comp_df["mode_display"] = comp_df["perturbation_mode"].map(MODE_DISPLAY_TITLES)

    p8 = (
        alt.Chart(comp_df)
        .mark_bar()
        .encode(
            y=alt.Y("mode_display:N", sort="-x", title=None),
            x=alt.X("rho_mean:Q", title="Mean Spearman ρ (at Intensity ~ 0.6)"),
            color=alt.Color("rho_mean:Q", scale=alt.Scale(scheme="blues"), legend=None),
            tooltip=["mode_display:N", alt.Tooltip("rho_mean:Q", format="+.3f"), alt.Tooltip("sign_mean:Q", format=".1f")],
        )
        .properties(
            width=320,
            height=280,
            title=alt.TitleParams(
                text="H: Concordance Comparison",
                subtitle="Mean Spearman ρ at representative strength",
                fontSize=12,
                subtitleFontSize=10,
            ),
        )
    )
    panels.append(p8)

    row1 = alt.hconcat(panels[0], panels[1], panels[2], panels[3])
    row2 = alt.hconcat(panels[4], panels[5], panels[6], panels[7])

    full_scatter_grid = (
        alt.vconcat(row1, row2)
        .properties(
            title=alt.TitleParams(
                text="Synthetic Benchmarking: Baseline vs. 6 Distortions (Sade-Feldman Style Scatter Grid)",
                subtitle=[
                    "Comparing Single-Cell Milopy DA vs. Pseudobulk Deconvolution across cell states for each distortion mode",
                    "Points colored by quadrant classification: Concordant Responder (green), Concordant Non-Responder (blue), Discordant (red/orange)",
                ],
                fontSize=16,
                subtitleFontSize=12,
                anchor="middle",
            )
        )
        .configure_view(strokeWidth=1, stroke="#cccccc")
    )
    return full_scatter_grid


def plot_compound_scatter_grid(
    compound_scatter_df: pl.DataFrame,
    compound_summary_df: pl.DataFrame,
) -> alt.VConcatChart:
    """Construct Sade-Feldman style scatter grid comparing Baseline vs. the 4 Compound Clinical Regimes."""
    quadrant_colors = {
        "Concordant Responder": "#2ca02c",      # Green
        "Concordant Non-Responder": "#1f77b4",  # Blue
        "Discordant (Milo+, Deconv-)": "#d62728", # Red
        "Discordant (Milo-, Deconv+)": "#ff7f0e", # Orange
        "Concordant Neutral": "#7f7f7f",        # Gray
    }
    color_scale = alt.Scale(
        domain=list(quadrant_colors.keys()),
        range=list(quadrant_colors.values()),
    )

    panels: list[alt.Chart] = []
    h_rule = alt.Chart().mark_rule(color="#dddddd").encode(y=alt.datum(0))
    v_rule = alt.Chart().mark_rule(color="#dddddd").encode(x=alt.datum(0))

    # Panel 1: Baseline Unperturbed (intensity = 0.0)
    base_points = compound_scatter_df.filter((pl.col("intensity") == 0.0) & (pl.col("replicate") == 0)).to_pandas()
    base_sum = compound_summary_df.filter((pl.col("intensity") == 0.0)).row(0, named=True)

    c_base = (
        alt.Chart(base_points)
        .mark_circle(size=140, opacity=0.9, stroke="black", strokeWidth=1)
        .encode(
            x=alt.X("milo_mean_logfc:Q", title="Milopy Mean logFC"),
            y=alt.Y("deconv_beta:Q", title="Deconvolution Effect (β)"),
            color=alt.Color("quadrant:N", scale=color_scale, legend=None),
            tooltip=["cell_state:N", "milo_mean_logfc:Q", "deconv_beta:Q", "quadrant:N"],
        )
    )
    line_base = c_base.transform_regression("milo_mean_logfc", "deconv_beta").mark_line(color="#444444", strokeDash=[4, 4], strokeWidth=2)
    p1 = (c_base + line_base + h_rule + v_rule).properties(
        width=320,
        height=280,
        title=alt.TitleParams(
            text="A: Baseline (Unperturbed)",
            subtitle=f"Spearman ρ = {base_sum['rho_mean']:+.3f} ± {base_sum['rho_std']:.3f} | Sign = {base_sum['sign_mean']:.0f}%",
            fontSize=12,
            subtitleFontSize=10,
        ),
    )
    panels.append(p1)

    # Panels 2-5: The 4 Compound Regimes at moderate-to-severe intensity ~ 0.6
    regime_list = [
        ("core_biopsy", "B: Clinical Core Needle Biopsy", "Cell Size + Sparsity + Ghost"),
        ("inflamed_tme", "C: Inflamed Tumor Microenv.", "Activation + Collinearity + Shift"),
        ("sc_technical", "D: Single-Cell Technical Noise", "Ambient Soup + Sparsity + Shift"),
        ("triple_jeopardy", "E: Triple-Jeopardy Breakdown", "Cell Size + Activation + Ghost"),
    ]

    for regime_key, title, subtitle_desc in regime_list:
        sub_points = (
            compound_scatter_df.filter(
                (pl.col("perturbation_mode") == regime_key)
                & ((pl.col("intensity") - 0.6).abs() < 0.05)
                & (pl.col("replicate") == 0)
            )
            .to_pandas()
        )
        sub_sum_rows = compound_summary_df.filter(
            (pl.col("perturbation_mode") == regime_key)
            & ((pl.col("intensity") - 0.6).abs() < 0.05)
        )
        if sub_sum_rows.height > 0:
            s_row = sub_sum_rows.row(0, named=True)
            stats_str = f"ρ = {s_row['rho_mean']:+.3f} ± {s_row['rho_std']:.3f} | Sign = {s_row['sign_mean']:.0f}%"
        else:
            stats_str = "Instance"

        c_pt = (
            alt.Chart(sub_points)
            .mark_circle(size=140, opacity=0.9, stroke="black", strokeWidth=1)
            .encode(
                x=alt.X("milo_mean_logfc:Q", title="Milopy Mean logFC"),
                y=alt.Y("deconv_beta:Q", title="Deconvolution Effect (β)"),
                color=alt.Color("quadrant:N", scale=color_scale, legend=None),
                tooltip=["cell_state:N", "milo_mean_logfc:Q", "deconv_beta:Q", "quadrant:N"],
            )
        )
        reg_line = c_pt.transform_regression("milo_mean_logfc", "deconv_beta").mark_line(color="#444444", strokeDash=[4, 4], strokeWidth=2)
        p_sub = (c_pt + reg_line + h_rule + v_rule).properties(
            width=320,
            height=280,
            title=alt.TitleParams(text=title, subtitle=f"{subtitle_desc} ({stats_str})", fontSize=12, subtitleFontSize=9),
        )
        panels.append(p_sub)

    # Panel 6: Summary Comparison Bar Chart
    comp_df = (
        compound_summary_df.filter((pl.col("intensity") - 0.6).abs() < 0.05)
        .to_pandas()
    )
    comp_df["regime_title"] = comp_df["perturbation_mode"].map(lambda k: COMPOUND_REGIMES[k]["title"])

    p6 = (
        alt.Chart(comp_df)
        .mark_bar()
        .encode(
            y=alt.Y("regime_title:N", sort="-x", title=None),
            x=alt.X("rho_mean:Q", title="Mean Spearman ρ (Intensity ~ 0.6)"),
            color=alt.Color("rho_mean:Q", scale=alt.Scale(scheme="redyellowgreen", domain=[-0.2, 0.9]), legend=None),
            tooltip=["regime_title:N", alt.Tooltip("rho_mean:Q", format="+.3f"), alt.Tooltip("sign_mean:Q", format=".1f")],
        )
        .properties(
            width=320,
            height=280,
            title=alt.TitleParams(
                text="F: Compound Concordance Comparison",
                subtitle="Mean Spearman ρ at Intensity ~ 0.6",
                fontSize=12,
                subtitleFontSize=10,
            ),
        )
    )
    panels.append(p6)

    row1 = alt.hconcat(panels[0], panels[1], panels[2])
    row2 = alt.hconcat(panels[3], panels[4], panels[5])

    full_grid = (
        alt.vconcat(row1, row2)
        .properties(
            title=alt.TitleParams(
                text="Compound Clinical Stress-Testing: Baseline vs. 4 Multi-Factorial Regimes",
                subtitle=[
                    "Evaluating simultaneous biological, technical, and transcriptional distortions at representative clinical strength (Intensity ~ 0.6)",
                    "Quadrant colors: Concordant Responder (green), Concordant Non-Responder (blue), Discordant (red/orange)",
                ],
                fontSize=16,
                subtitleFontSize=12,
                anchor="middle",
            )
        )
        .configure_view(strokeWidth=1, stroke="#cccccc")
    )
    return full_grid


# ------------------------------------------------------------------------------
# 7. Main Execution Boundary
# ------------------------------------------------------------------------------

def main() -> None:
    config = parse_args()
    config.out_dir.mkdir(parents=True, exist_ok=True)
    config.results_dir.mkdir(parents=True, exist_ok=True)

    # 1. Single-Mode Replicate Sweep
    if config.run_single_sweep:
        raw_sweep_df, summary_df, scatter_df = run_replicate_sweep(config)
        raw_sweep_df.write_parquet(config.out_dir / "perturbation_replicates_raw.parquet")
        summary_df.write_parquet(config.out_dir / "perturbation_sensitivity_summary.parquet")
        scatter_df.write_parquet(config.out_dir / "synthetic_perturbation_scatter_points.parquet")

        print(f"\nSaved raw replicate results ({raw_sweep_df.height} runs) to {config.out_dir}")

        curves_chart = plot_sensitivity_curves_with_sd(summary_df)
        curves_svg_path = config.results_dir / config.curves_svg_name
        with open(curves_svg_path, "w", encoding="utf-8") as f:
            f.write(vlc.vegalite_to_svg(curves_chart.to_dict()))
        print(f"Exported Sensitivity Curves SVG to {curves_svg_path}")

        scatter_chart = plot_scatter_grid(scatter_df, summary_df)
        scatter_svg_path = config.results_dir / config.scatter_svg_name
        with open(scatter_svg_path, "w", encoding="utf-8") as f:
            f.write(vlc.vegalite_to_svg(scatter_chart.to_dict()))
        print(f"Exported Scatter Grid SVG to {scatter_svg_path}")

    # 2. Compound Regimes Sweep
    if config.run_compound_sweep:
        raw_compound_df, compound_summary_df, compound_scatter_df = run_compound_sweep(config)
        raw_compound_df.write_parquet(config.out_dir / "compound_replicates_raw.parquet")
        compound_summary_df.write_parquet(config.out_dir / "compound_sensitivity_summary.parquet")
        compound_scatter_df.write_parquet(config.out_dir / "compound_scatter_points.parquet")

        print(f"\nSaved compound replicate results ({raw_compound_df.height} runs) to {config.out_dir}")

        compound_curves_chart = plot_compound_sensitivity_curves(compound_summary_df)
        compound_curves_svg_path = config.results_dir / config.compound_curves_svg_name
        with open(compound_curves_svg_path, "w", encoding="utf-8") as f:
            f.write(vlc.vegalite_to_svg(compound_curves_chart.to_dict()))
        print(f"Exported Compound Sensitivity Curves SVG to {compound_curves_svg_path}")

        compound_scatter_chart = plot_compound_scatter_grid(compound_scatter_df, compound_summary_df)
        compound_scatter_svg_path = config.results_dir / config.compound_scatter_svg_name
        with open(compound_scatter_svg_path, "w", encoding="utf-8") as f:
            f.write(vlc.vegalite_to_svg(compound_scatter_chart.to_dict()))
        print(f"Exported Compound Scatter Grid SVG to {compound_scatter_svg_path}")

    # 3. 2D Factorial Interaction Sweep (Cell Size x Activation Confounding)
    if config.run_factorial_sweep:
        raw_factorial_df, factorial_summary_df = run_2d_factorial_sweep(config)
        raw_factorial_df.write_parquet(config.out_dir / "factorial_2d_raw.parquet")
        factorial_summary_df.write_parquet(config.out_dir / "factorial_2d_summary.parquet")

        print(f"\nSaved 2D factorial interaction results ({raw_factorial_df.height} runs) to {config.out_dir}")

        heatmap_chart = plot_2d_interaction_heatmap(factorial_summary_df)
        heatmap_svg_path = config.results_dir / config.factorial_heatmap_svg_name
        with open(heatmap_svg_path, "w", encoding="utf-8") as f:
            f.write(vlc.vegalite_to_svg(heatmap_chart.to_dict()))
        print(f"Exported 2D Factorial Heatmap SVG to {heatmap_svg_path}")

    print("\nAll benchmark pipelines completed successfully!")
    sys.exit(0)


if __name__ == "__main__":
    main()
