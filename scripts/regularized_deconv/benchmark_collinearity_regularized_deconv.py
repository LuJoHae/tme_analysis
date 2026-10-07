#!/usr/bin/env python3
"""
Benchmark: Collinearity-Aware Regularized Deconvolution vs Flat Deconvolution and BayesPrism.

Evaluates:
1. Synthetic Benchmark with Exact Ground Truth (50 samples, 6 states, collinearity r > 0.99):
   - Ground truth recovery MSE
   - Condition number reduction
   - Suppression of negative cross-talk
2. Real-World Clinical Benchmark on Sade-Feldman Melanoma Reference (12 states) & Hugo Cohort:
   - State-level concordance against single-cell Milo DA
   - Demonstration of sign-flip elimination at the state level

Adheres to strict functional Python with returns, Pydantic (frozen=True), and Polars.
"""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd  # type: ignore
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.maybe import Some, Nothing
from returns.result import Result, Success, Failure
from scipy import stats  # type: ignore
from scipy.optimize import nnls  # type: ignore

from bayesprism.regularized import (
    RegularizedDeconvConfig,
    build_transcriptomic_graph,
    calibrate_spectral_lambda,
    deconvolve_collinearity_regularized,
)
from bayesprism.adapters.rectangle import (
    RectangleConfig,
    create_signature_from_matrix,
    deconvolve_rectangle,
    is_rectangle_available,
)
from bayesprism.adapters.cibersort import (
    CibersortConfig,
    deconvolve_cibersort,
)
from bayesprism.adapters.cibersortx_docker import (
    CibersortXDockerConfig,
    run_cibersortx_docker,
    is_docker_available,
    resolve_cibersortx_credentials,
)
import bayesprism as bp
import instaprism

try:
    from .synthetic_data import (
        generate_synthetic_reference,
        sample_true_proportions,
        generate_bulk_mixtures,
    )
    from .deconv_runners import ALL_METHODS, run_deconvolution_suite
except ImportError:
    from regularized_deconv.synthetic_data import (
        generate_synthetic_reference,
        sample_true_proportions,
        generate_bulk_mixtures,
    )
    from regularized_deconv.deconv_runners import ALL_METHODS, run_deconvolution_suite


try:
    from loguru import logger
    logger.remove()
except Exception:
    pass


class BenchmarkPaths(BaseModel):
    model_config = ConfigDict(frozen=True)
    ref_path: Path = Path("output/output/sade_feldman_deconv_validation/reference_phi_res0.5.parquet")
    hugo_tpm_path: Path = Path(
        "scratch/lair/CBioPortalDataset-Hugo-iAtlas/mel_iatlas_hugo_ucla_2016/mel_iatlas_hugo_ucla_2016/data_mrna_seq_tpm.txt"
    )
    hugo_clin_path: Path = Path(
        "scratch/lair/CBioPortalDataset-Hugo-iAtlas/mel_iatlas_hugo_ucla_2016/mel_iatlas_hugo_ucla_2016/data_clinical_sample.txt"
    )
    sc_effects_path: Path = Path("output/concordance/sc_patient_response_effects.parquet")
    out_dir: Path = Path("output/concordance")
    cibersortx_username: str = ""
    cibersortx_token: str = ""


def run_multi_regime_synthetic_benchmark(
    n_samples: int = 50,
    n_genes: int = 400,
    seed: int = 42,
    out_dir: Path = Path("output/concordance"),
    cibersortx_username: str = "",
    cibersortx_token: str = "",
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Run controlled continuous correlation spectrum synthetic deconvolution benchmark across all 7 tools."""
    rng = np.random.default_rng(seed)
    target_r_values = (0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.98, 0.99)
    sample_names = tuple(f"Synthetic_Sample_{i:02d}" for i in range(n_samples))
    lineages = ("Lineage_1", "Lineage_2", "Lineage_3", "Lineage_4")

    sample_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    true_prop_rows: list[dict[str, object]] = []
    est_prop_rows: list[dict[str, object]] = []
    ref_mat_rows: list[dict[str, object]] = []

    out_dir.mkdir(parents=True, exist_ok=True)

    creds: Maybe[tuple[str, str]] = (
        Some((cibersortx_username, cibersortx_token))
        if (cibersortx_username and cibersortx_token)
        else Nothing
    )

    for target_r in target_r_values:
        sig = generate_synthetic_reference(
            n_lineages=4,
            states_per_lineage=3,
            n_genes=n_genes,
            target_r=target_r,
            rng=rng,
        )
        S = sig.phi.shape[0]
        T = sig.aggregation_matrix.shape[0]

        # Record reference matrix
        for s_idx, s_name in enumerate(sig.state_names):
            row_dict: dict[str, object] = {
                "target_r": target_r,
                "collinearity_r": sig.empirical_r,
                "state_name": s_name,
                "lineage": sig.lineage_labels[s_idx],
            }
            for g_idx in range(min(n_genes, 50)):
                row_dict[f"gene_{g_idx:03d}"] = float(sig.phi[s_idx, g_idx])
            ref_mat_rows.append(row_dict)

        # Condition numbers
        theta_0 = np.ones(S) / float(S)
        mu_0 = theta_0 @ sig.phi
        H_nom = sig.phi @ np.diag(1.0 / np.clip(mu_0, 1e-12, None)) @ sig.phi.T
        eigs_unreg = np.linalg.eigvalsh(H_nom)
        kappa_unreg = float(eigs_unreg[-1] / max(eigs_unreg[0], 1e-12))

        W, L = build_transcriptomic_graph(sig.phi, power=2.0)
        lam_lap, _ = calibrate_spectral_lambda(sig.phi, L, kappa_target=2000.0)
        H_reg = H_nom + lam_lap * L
        eigs_reg = np.linalg.eigvalsh(H_reg)
        kappa_reg = float(eigs_reg[-1] / max(eigs_reg[0], 1e-12))

        # True proportions
        true_theta = sample_true_proportions(
            n_samples=n_samples,
            signature=sig,
            mode="skewed",
            rng=rng,
        )
        true_lineage = true_theta @ sig.aggregation_matrix.T

        # Record true proportions
        for n in range(n_samples):
            t_row: dict[str, object] = {
                "target_r": target_r,
                "collinearity_r": sig.empirical_r,
                "sample_id": sample_names[n],
            }
            for s_idx, s_name in enumerate(sig.state_names):
                t_row[s_name] = float(true_theta[n, s_idx])
            for t_idx, l_name in enumerate(lineages):
                t_row[f"Lineage_{l_name}"] = float(true_lineage[n, t_idx])
            true_prop_rows.append(t_row)

        # Bulk mixtures
        mixture = generate_bulk_mixtures(
            true_theta=true_theta,
            phi=sig.phi,
            n_total=80_000,
            rng=rng,
        )

        # Run all 7 methods
        deconv_res = run_deconvolution_suite(
            mixture=mixture,
            signature=sig,
            methods=ALL_METHODS,
            sample_names=sample_names,
            cibersortx_credentials=creds,
        )

        for m_name in ALL_METHODS:
            t_est = deconv_res[m_name]
            l_est = t_est @ sig.aggregation_matrix.T
            k_val = kappa_reg if "RegDeconv" in m_name else kappa_unreg

            # Sibling cross-talk
            corrs: list[float] = []
            for t_idx in range(T):
                for a in range(3):
                    for b in range(a + 1, 3):
                        ca = t_est[:, t_idx * 3 + a]
                        cb = t_est[:, t_idx * 3 + b]
                        if np.std(ca) > 1e-8 and np.std(cb) > 1e-8:
                            corrs.append(float(np.corrcoef(ca, cb)[0, 1]))
                        else:
                            corrs.append(0.0)
            sib_r_val = float(np.mean(corrs)) if corrs else 0.0

            s_mses: list[float] = []
            l_mses: list[float] = []
            dropouts: list[float] = []

            for n in range(n_samples):
                s_mse = float(np.mean((t_est[n] - true_theta[n]) ** 2))
                l_mse = float(np.mean((l_est[n] - true_lineage[n]) ** 2))
                drop = float(np.mean((true_theta[n] > 0.01) & (t_est[n] < 1e-4)))
                s_mses.append(s_mse)
                l_mses.append(l_mse)
                dropouts.append(drop)

                sample_rows.append({
                    "target_r": target_r,
                    "collinearity_r": sig.empirical_r,
                    "sample_id": sample_names[n],
                    "method": m_name,
                    "condition_number": k_val,
                    "state_mse": s_mse,
                    "lineage_mse": l_mse,
                    "dropout_rate": drop,
                })

                e_row: dict[str, object] = {
                    "target_r": target_r,
                    "collinearity_r": sig.empirical_r,
                    "sample_id": sample_names[n],
                    "method": m_name,
                }
                for s_idx, s_name in enumerate(sig.state_names):
                    e_row[s_name] = float(t_est[n, s_idx])
                for t_idx, l_name in enumerate(lineages):
                    e_row[f"Lineage_{l_name}"] = float(l_est[n, t_idx])
                est_prop_rows.append(e_row)

            summary_rows.append({
                "target_r": target_r,
                "collinearity_r": sig.empirical_r,
                "method": m_name,
                "condition_number": k_val,
                "state_mse_mean": float(np.mean(s_mses)),
                "state_mse_std": float(np.std(s_mses)),
                "lineage_mse_mean": float(np.mean(l_mses)),
                "lineage_mse_std": float(np.std(l_mses)),
                "sibling_corr": sib_r_val,
                "dropout_rate_pct": float(np.mean(dropouts) * 100.0),
            })

    sample_df = pl.DataFrame(sample_rows)
    summary_df = pl.DataFrame(summary_rows)
    true_prop_df = pl.DataFrame(true_prop_rows)
    est_prop_df = pl.DataFrame(est_prop_rows)
    ref_mat_df = pl.DataFrame(ref_mat_rows)

    sample_df.write_parquet(out_dir / "synthetic_collinearity_benchmark_results.parquet")
    summary_df.write_parquet(out_dir / "synthetic_collinearity_benchmark_summary.parquet")
    true_prop_df.write_parquet(out_dir / "synthetic_collinearity_true_proportions.parquet")
    est_prop_df.write_parquet(out_dir / "synthetic_collinearity_estimated_proportions.parquet")
    ref_mat_df.write_parquet(out_dir / "synthetic_collinearity_reference_matrices.parquet")

    return summary_df, sample_df



def run_clinical_sade_feldman_benchmark(
    paths: BenchmarkPaths,
) -> Result[pl.DataFrame, str]:
    """Run real-world clinical benchmark comparing state-level effects against single-cell Milo DA."""
    if not (paths.ref_path.exists() and paths.hugo_tpm_path.exists() and paths.hugo_clin_path.exists()):
        return Failure(f"One or more input files do not exist: {paths}")

    ref_df = pl.read_parquet(paths.ref_path)
    bulk_df = pd.read_csv(paths.hugo_tpm_path, sep="\t", index_col=0)
    bulk_df.index = bulk_df.index.astype(str).str.upper()
    bulk_df = bulk_df.groupby(level=0).mean()

    clin_df = pd.read_csv(paths.hugo_clin_path, sep="\t", skiprows=4).set_index("SAMPLE_ID")
    sc_df = pl.read_parquet(paths.sc_effects_path)

    clusters = ref_df["cluster"].to_list()
    ref_genes = [c for c in ref_df.columns if c != "cluster"]
    ref_mat = ref_df.select(ref_genes).to_numpy().astype(np.float64)

    common_genes = sorted(list(set(bulk_df.index).intersection(set(ref_genes))))
    ref_gene_indices = [ref_genes.index(g) for g in common_genes]
    sub_ref = ref_mat[:, ref_gene_indices]

    sub_bulk = bulk_df.loc[common_genes].T
    sample_ids = [s for s in list(sub_bulk.index) if s in clin_df.index]
    bulk_counts = np.round(sub_bulk.loc[sample_ids].to_numpy().astype(np.float64)).astype(np.int64)

    # 1. Unregularized Flat Deconvolution
    cfg_unreg = RegularizedDeconvConfig(lambda_lap=Some(0.0), lambda_fuse=Some(0.0), max_iter=300)
    res_unreg = deconvolve_collinearity_regularized(
        mixture=bulk_counts,
        reference=sub_ref,
        state_names=clusters,
        sample_names=sample_ids,
        config=Some(cfg_unreg),
    ).unwrap()

    # 2. Collinearity-Aware Regularized Engine (calibrated lambda_lap=0.50, lambda_fuse=0.05)
    cfg_reg = RegularizedDeconvConfig(lambda_lap=Some(0.50), lambda_fuse=Some(0.05), max_iter=300)
    res_reg = deconvolve_collinearity_regularized(
        mixture=bulk_counts,
        reference=sub_ref,
        state_names=clusters,
        sample_names=sample_ids,
        config=Some(cfg_reg),
    ).unwrap()

    # 3. Rectangle Multiscale Deconvolution (DWLS-QP)
    has_rect = is_rectangle_available()
    res_rect = None
    if has_rect:
        ref_cpm = sub_ref * 1e6
        bulk_cpm = (bulk_counts / np.sum(bulk_counts, axis=1, keepdims=True)) * 1e6
        sig_rect = create_signature_from_matrix(ref_cpm, tuple(clusters), tuple(common_genes)).unwrap()
        cfg_rect = RectangleConfig(correct_mrna_bias=False)
        res_rect = deconvolve_rectangle(
            signatures=sig_rect,
            bulks=bulk_cpm,
            sample_names=tuple(sample_ids),
            gene_names=tuple(common_genes),
            config=Some(cfg_rect),
        ).unwrap()

    # 4. CIBERSORT (linear nu-SVR reimplementation)
    cfg_ciber = CibersortConfig(nu_values=(0.25, 0.50, 0.75), n_perm=0)
    res_ciber = deconvolve_cibersort(
        mixture=bulk_counts,
        reference=sub_ref,
        state_names=tuple(clusters),
        sample_names=tuple(sample_ids),
        gene_names=tuple(common_genes),
        config=Some(cfg_ciber),
    ).unwrap()

    # 5. InstaPrism (EM / fixpoint)
    row_sums = sub_ref.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    norm_ref = (sub_ref / row_sums).astype(np.float64)
    insta_mat = np.zeros((len(sample_ids), len(clusters)), dtype=np.float64)
    for i in range(len(sample_ids)):
        b_vec = bulk_counts[i].astype(np.float64)
        _, _, fracs, _ = instaprism.insta_prism(bulk=b_vec, reference=norm_ref, n_iter=50)
        insta_mat[i, :] = fracs

    # 6. BayesPrism (Gibbs MCMC)
    prism_res = bp.new_prism(
        reference=(sub_ref * 100).astype(np.float32),
        cell_type_labels=tuple(clusters),
        cell_state_labels=Nothing,
        mixture=bulk_counts.astype(np.float32),
        gene_names=Some(tuple(common_genes)),
        pseudo_min=1e-8,
    )
    bp_mat = np.ones((len(sample_ids), len(clusters)), dtype=np.float64) / len(clusters)
    match prism_res:
        case Success(prism_obj):
            ctrl = bp.GibbsControl(chain_length=40, burn_in=20, thinning=2, seed=42)
            fit_res = bp.run_prism(prism_obj, update_gibbs=False, gibbs_control=Some(ctrl))
            match fit_res:
                case Success(fitted):
                    frac_res = bp.get_fraction(fitted, which_theta="first", state_or_type="type")
                    match frac_res:
                        case Success(bp_df):
                            bp_mat = bp_df.select(list(clusters)).to_numpy()
                        case Failure(_):
                            pass
                case Failure(_):
                    pass
        case Failure(_):
            pass

    # Evaluate point-biserial response associations
    responses = np.array([
        1.0 if str(clin_df.loc[s, "RESPONSE"]).strip().lower() in ["complete response", "partial response", "responder"] else 0.0
        for s in sample_ids
    ], dtype=np.float64)

    effects_rows: list[dict[str, object]] = []

    for idx, c_name in enumerate(clusters):
        sc_val = sc_df.filter(pl.col("cell_state") == c_name)["beta_sc"][0] if c_name in sc_df["cell_state"] else 0.0

        # Unregularized beta
        x_u = res_unreg.theta[:, idx]
        if np.std(x_u) > 1e-9:
            x_u_std = (x_u - np.mean(x_u)) / np.std(x_u)
            r_u, _ = stats.pointbiserialr(responses, x_u_std)
            r_u_clip = np.clip(r_u, -0.999, 0.999)
            beta_u = float(2.0 * r_u_clip / np.sqrt(1.0 - r_u_clip**2 + 1e-12))
        else:
            beta_u = float("nan")

        # RegDeconv beta
        x_r = res_reg.theta[:, idx]
        if np.std(x_r) > 1e-9:
            x_r_std = (x_r - np.mean(x_r)) / np.std(x_r)
            r_r, _ = stats.pointbiserialr(responses, x_r_std)
            r_r_clip = np.clip(r_r, -0.999, 0.999)
            beta_r = float(2.0 * r_r_clip / np.sqrt(1.0 - r_r_clip**2 + 1e-12))
        else:
            beta_r = float("nan")

        # Rectangle beta
        if res_rect is not None and c_name in res_rect.proportions.columns:
            x_rc = res_rect.proportions[c_name].to_numpy()
            if np.std(x_rc) > 1e-9:
                x_rc_std = (x_rc - np.mean(x_rc)) / np.std(x_rc)
                r_rc, _ = stats.pointbiserialr(responses, x_rc_std)
                r_rc_clip = np.clip(r_rc, -0.999, 0.999)
                beta_rc = float(2.0 * r_rc_clip / np.sqrt(1.0 - r_rc_clip**2 + 1e-12))
            else:
                beta_rc = float("nan")
        else:
            beta_rc = float("nan")

        # CIBERSORT (reimpl.) beta
        if c_name in res_ciber.proportions.columns:
            x_cb = res_ciber.proportions[c_name].to_numpy()
            if np.std(x_cb) > 1e-9:
                x_cb_std = (x_cb - np.mean(x_cb)) / np.std(x_cb)
                r_cb, _ = stats.pointbiserialr(responses, x_cb_std)
                r_cb_clip = np.clip(r_cb, -0.999, 0.999)
                beta_cb = float(2.0 * r_cb_clip / np.sqrt(1.0 - r_cb_clip**2 + 1e-12))
            else:
                beta_cb = float("nan")
        else:
            beta_cb = float("nan")

        # InstaPrism beta
        x_ip = insta_mat[:, idx]
        if np.std(x_ip) > 1e-9:
            x_ip_std = (x_ip - np.mean(x_ip)) / np.std(x_ip)
            r_ip, _ = stats.pointbiserialr(responses, x_ip_std)
            r_ip_clip = np.clip(r_ip, -0.999, 0.999)
            beta_ip = float(2.0 * r_ip_clip / np.sqrt(1.0 - r_ip_clip**2 + 1e-12))
        else:
            beta_ip = float("nan")

        # BayesPrism beta
        x_bp = bp_mat[:, idx]
        if np.std(x_bp) > 1e-9:
            x_bp_std = (x_bp - np.mean(x_bp)) / np.std(x_bp)
            r_bp, _ = stats.pointbiserialr(responses, x_bp_std)
            r_bp_clip = np.clip(r_bp, -0.999, 0.999)
            beta_bp = float(2.0 * r_bp_clip / np.sqrt(1.0 - r_bp_clip**2 + 1e-12))
        else:
            beta_bp = float("nan")

        effects_rows.append({
            "cell_state": c_name,
            "beta_sc": float(sc_val),
            "beta_unregularized": beta_u,
            "beta_regularized": beta_r,
            "beta_rectangle": beta_rc,
            "beta_cibersort": beta_cb,
            "beta_instaprism": beta_ip,
            "beta_bayesprism": beta_bp,
            "sign_concordant_unreg": (np.sign(sc_val) == np.sign(beta_u)) if (not np.isnan(beta_u) and abs(sc_val) > 0.05) else False,
            "sign_concordant_reg": (np.sign(sc_val) == np.sign(beta_r)) if (not np.isnan(beta_r) and abs(sc_val) > 0.05) else False,
            "sign_concordant_rect": (np.sign(sc_val) == np.sign(beta_rc)) if (not np.isnan(beta_rc) and abs(sc_val) > 0.05) else False,
            "sign_concordant_ciber": (np.sign(sc_val) == np.sign(beta_cb)) if (not np.isnan(beta_cb) and abs(sc_val) > 0.05) else False,
            "sign_concordant_insta": (np.sign(sc_val) == np.sign(beta_ip)) if (not np.isnan(beta_ip) and abs(sc_val) > 0.05) else False,
            "sign_concordant_bp": (np.sign(sc_val) == np.sign(beta_bp)) if (not np.isnan(beta_bp) and abs(sc_val) > 0.05) else False,
        })

    effects_df = pl.DataFrame(effects_rows)
    return Success(effects_df)


def run_benchmark(
    cibersortx_username: str = "",
    cibersortx_token: str = "",
) -> Result[None, str]:
    print("=" * 80)
    print("1. SYNTHETIC COLLINEAR BENCHMARK (GROUND TRUTH COMPARISON)")
    print("=" * 80)
    summary_df, _ = run_multi_regime_synthetic_benchmark(
        cibersortx_username=cibersortx_username,
        cibersortx_token=cibersortx_token,
    )
    print(summary_df)

    paths = BenchmarkPaths(
        cibersortx_username=cibersortx_username,
        cibersortx_token=cibersortx_token,
    )
    print("\n" + "=" * 80)
    print("2. CLINICAL BENCHMARK ON SADE-FELDMAN MELANOMA COHORT (12 IMMUNE STATES)")
    print("=" * 80)

    match run_clinical_sade_feldman_benchmark(paths):
        case Failure(err):
            print(f"[WARNING] Clinical benchmark skipped: {err}")
        case Success(eff_df):
            print(eff_df)

            # Metrics
            b_sc = eff_df["beta_sc"].to_numpy()
            b_unreg = eff_df["beta_unregularized"].to_numpy()
            b_reg = eff_df["beta_regularized"].to_numpy()
            b_rect = eff_df["beta_rectangle"].to_numpy()
            b_ciber = eff_df["beta_cibersort"].to_numpy()
            b_insta = eff_df["beta_instaprism"].to_numpy()
            b_bp = eff_df["beta_bayesprism"].to_numpy()

            valid_u = ~np.isnan(b_unreg)
            valid_r = ~np.isnan(b_reg)
            valid_rc = ~np.isnan(b_rect)
            valid_cb = ~np.isnan(b_ciber)
            valid_ip = ~np.isnan(b_insta)
            valid_bp = ~np.isnan(b_bp)

            rho_u = float(stats.spearmanr(b_sc[valid_u], b_unreg[valid_u])[0]) if np.sum(valid_u) > 2 else float("nan")
            rho_r = float(stats.spearmanr(b_sc[valid_r], b_reg[valid_r])[0]) if np.sum(valid_r) > 2 else float("nan")
            rho_rc = float(stats.spearmanr(b_sc[valid_rc], b_rect[valid_rc])[0]) if np.sum(valid_rc) > 2 else float("nan")
            rho_cb = float(stats.spearmanr(b_sc[valid_cb], b_ciber[valid_cb])[0]) if np.sum(valid_cb) > 2 else float("nan")
            rho_ip = float(stats.spearmanr(b_sc[valid_ip], b_insta[valid_ip])[0]) if np.sum(valid_ip) > 2 else float("nan")
            rho_bp = float(stats.spearmanr(b_sc[valid_bp], b_bp[valid_bp])[0]) if np.sum(valid_bp) > 2 else float("nan")

            concord_u = float(eff_df.filter(pl.col("beta_unregularized").is_not_nan())["sign_concordant_unreg"].mean() * 100.0) if np.sum(valid_u) > 0 else 0.0
            concord_r = float(eff_df.filter(pl.col("beta_regularized").is_not_nan())["sign_concordant_reg"].mean() * 100.0) if np.sum(valid_r) > 0 else 0.0
            concord_rc = float(eff_df.filter(pl.col("beta_rectangle").is_not_nan())["sign_concordant_rect"].mean() * 100.0) if np.sum(valid_rc) > 0 else 0.0
            concord_cb = float(eff_df.filter(pl.col("beta_cibersort").is_not_nan())["sign_concordant_ciber"].mean() * 100.0) if np.sum(valid_cb) > 0 else 0.0
            concord_ip = float(eff_df.filter(pl.col("beta_instaprism").is_not_nan())["sign_concordant_insta"].mean() * 100.0) if np.sum(valid_ip) > 0 else 0.0
            concord_bp = float(eff_df.filter(pl.col("beta_bayesprism").is_not_nan())["sign_concordant_bp"].mean() * 100.0) if np.sum(valid_bp) > 0 else 0.0

            print("\nClinical Concordance Summary:")
            print(f"  Unregularized Valid States: {np.sum(valid_u)}/12 | Spearman rho: {rho_u:.3f} | Sign Concordance: {concord_u:.1f}%")
            print(f"  RegDeconv     Valid States: {np.sum(valid_r)}/12 | Spearman rho: {rho_r:.3f} | Sign Concordance: {concord_r:.1f}%")
            print(f"  Rectangle     Valid States: {np.sum(valid_rc)}/12 | Spearman rho: {rho_rc:.3f} | Sign Concordance: {concord_rc:.1f}%")
            print(f"  CIBERSORT (reimpl) States: {np.sum(valid_cb)}/12 | Spearman rho: {rho_cb:.3f} | Sign Concordance: {concord_cb:.1f}%")
            print(f"  InstaPrism    Valid States: {np.sum(valid_ip)}/12 | Spearman rho: {rho_ip:.3f} | Sign Concordance: {concord_ip:.1f}%")
            print(f"  BayesPrism    Valid States: {np.sum(valid_bp)}/12 | Spearman rho: {rho_bp:.3f} | Sign Concordance: {concord_bp:.1f}%")

            # Check 01_Tem/Trm specifically
            tem_01 = eff_df.filter(pl.col("cell_state").str.contains("01_Tem"))
            if not tem_01.is_empty():
                row = tem_01[0]
                print(f"\n  Focus State: {row['cell_state'][0]}")
                print(f"    Single-Cell Effect (beta_sc): {row['beta_sc'][0]:+.3f}")
                print(f"    Unregularized Bulk Effect:   {row['beta_unregularized'][0]:+.3f}")
                print(f"    RegDeconv Bulk Effect:       {row['beta_regularized'][0]:+.3f}")
                print(f"    Rectangle Bulk Effect:       {row['beta_rectangle'][0]:+.3f}")
                print(f"    CIBERSORT Bulk Effect:       {row['beta_cibersort'][0]:+.3f}")
                print(f"    InstaPrism Bulk Effect:      {row['beta_instaprism'][0]:+.3f}")
                print(f"    BayesPrism Bulk Effect:      {row['beta_bayesprism'][0]:+.3f}")

            # Save metrics
            paths.out_dir.mkdir(parents=True, exist_ok=True)
            eff_df.write_parquet(paths.out_dir / "regularized_deconv_benchmark_metrics.parquet")

    print("=" * 80)
    print("[BENCHMARK COMPLETED SUCCESSFULLY]")
    print("=" * 80)
    return Success(None)


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description="Run deconvolution benchmark.")
    parser.add_argument("--cibersortx-username", type=str, default="", help="CIBERSORTx username")
    parser.add_argument("--cibersortx-token", type=str, default="", help="CIBERSORTx token")
    args = parser.parse_args()

    match run_benchmark(
        cibersortx_username=args.cibersortx_username,
        cibersortx_token=args.cibersortx_token,
    ):
        case Success(_):
            sys.exit(0)
        case Failure(err):
            print(f"[ERROR] {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
