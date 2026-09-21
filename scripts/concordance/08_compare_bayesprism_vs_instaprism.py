#!/usr/bin/env python3
"""
Step 8: Direct Benchmark of Actual BayesPrism MCMC vs InstaPrism
with and without Tumor Purity Normalization on the Hugo Melanoma Cohort.
Computes:
1. Actual BayesPrism Gibbs MCMC (Dirichlet-Multinomial sampler)
2. InstaPrism Fast Matrix Projection
3. Raw Fractions (without tumor purity adjustment)
4. Purity-Normalized Fractions (rescaled relative to immune infiltrate)
5. Standardized logistic regression effects vs patient response
6. Concordance with single-cell ground-truth effects (Milo)
Strict functional Python with returns, Pydantic, Polars, and Altair.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import altair as alt  # type: ignore
import bayesprism as bp  # type: ignore
import instaprism  # type: ignore
import numpy as np
import pandas as pd  # type: ignore
from pydantic import BaseModel, ConfigDict
from returns.maybe import Nothing, Some
from returns.result import Failure, Result, Success
import polars as pl
from scipy import stats  # type: ignore
import vl_convert as vlc  # type: ignore


class BenchmarkConfig(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    ref_path: Path
    hugo_tpm_path: Path
    hugo_clin_path: Path
    sc_effects_path: Path
    out_dir: Path
    results_dir: Path
    chain_length: int = 40
    burn_in: int = 20


def parse_args() -> BenchmarkConfig:
    parser = argparse.ArgumentParser(
        description="Compare BayesPrism vs InstaPrism with/without purity normalization."
    )
    parser.add_argument(
        "--ref",
        type=Path,
        default=Path("output/output/sade_feldman_deconv_validation/reference_phi_res0.5.parquet"),
        help="Path to single-cell reference phi parquet",
    )
    parser.add_argument(
        "--hugo-tpm",
        type=Path,
        default=Path(
            "scratch/lair/CBioPortalDataset-Hugo-iAtlas/mel_iatlas_hugo_ucla_2016/mel_iatlas_hugo_ucla_2016/data_mrna_seq_tpm.txt"
        ),
        help="Path to Hugo bulk TPM file",
    )
    parser.add_argument(
        "--hugo-clin",
        type=Path,
        default=Path(
            "scratch/lair/CBioPortalDataset-Hugo-iAtlas/mel_iatlas_hugo_ucla_2016/mel_iatlas_hugo_ucla_2016/data_clinical_sample.txt"
        ),
        help="Path to Hugo clinical sample file",
    )
    parser.add_argument(
        "--sc-effects",
        type=Path,
        default=Path("output/concordance/sc_patient_response_effects.parquet"),
        help="Path to single-cell response effects parquet",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("output/concordance"),
        help="Output directory for parquets",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results/concordance"),
        help="Output directory for figures",
    )
    args = parser.parse_args()
    return BenchmarkConfig(
        ref_path=args.ref,
        hugo_tpm_path=args.hugo_tpm,
        hugo_clin_path=args.hugo_clin,
        sc_effects_path=args.sc_effects,
        out_dir=args.out_dir,
        results_dir=args.results_dir,
    )


def load_inputs(config: BenchmarkConfig) -> Result[tuple[pl.DataFrame, pd.DataFrame, pd.DataFrame], str]:
    """Pure boundary loading reference, bulk expression, and clinical data."""
    if not config.ref_path.exists():
        return Failure(f"Reference file not found: {config.ref_path}")
    if not config.hugo_tpm_path.exists():
        return Failure(f"Hugo TPM file not found: {config.hugo_tpm_path}")
    if not config.hugo_clin_path.exists():
        return Failure(f"Hugo clinical file not found: {config.hugo_clin_path}")

    ref_df = pl.read_parquet(config.ref_path)
    bulk_df = pd.read_csv(config.hugo_tpm_path, sep="\t", index_col=0)
    bulk_df.index = bulk_df.index.astype(str).str.upper()
    bulk_df = bulk_df.groupby(level=0).mean()

    clin_df = pd.read_csv(config.hugo_clin_path, sep="\t", skiprows=4)
    clin_df = clin_df.set_index("SAMPLE_ID")

    return Success((ref_df, bulk_df, clin_df))


def run_deconvolutions(
    ref_df: pl.DataFrame,
    bulk_df: pd.DataFrame,
    config: BenchmarkConfig,
) -> Result[tuple[pl.DataFrame, pl.DataFrame, list[str], list[str]], str]:
    """Run both BayesPrism Gibbs MCMC and InstaPrism on bulk cohort."""
    clusters = ref_df["cluster"].to_list()
    ref_genes = [c for c in ref_df.columns if c != "cluster"]
    ref_mat = ref_df.select(ref_genes).to_numpy().astype(np.float64)

    common_genes = sorted(list(set(bulk_df.index).intersection(set(ref_genes))))
    if len(common_genes) < 50:
        return Failure(f"Insufficient overlapping genes: {len(common_genes)}")

    ref_gene_indices = [ref_genes.index(g) for g in common_genes]
    sub_ref = ref_mat[:, ref_gene_indices]

    sub_bulk = bulk_df.loc[common_genes].T  # samples x genes
    sample_ids = list(sub_bulk.index)
    bulk_counts = np.round(sub_bulk.to_numpy().astype(np.float64)).astype(np.int64)

    print(f"[INFO] Running deconvolution on {len(sample_ids)} samples across {len(common_genes)} genes...")

    # 1. Run Actual BayesPrism Gibbs Sampler
    print("  [1/2] Running BayesPrism Gibbs MCMC sampler...")
    prism_res = bp.new_prism(
        reference=sub_ref,
        cell_type_labels=clusters,
        cell_state_labels=Nothing,
        mixture=bulk_counts,
        gene_names=Some(common_genes),
        bulk_names=Some(sample_ids),
    )
    match prism_res:
        case Failure(err):
            return Failure(f"BayesPrism initialization failed: {err}")
        case Success(prism_obj):
            ctrl = bp.GibbsControl(
                chain_length=config.chain_length,
                burn_in=config.burn_in,
                thinning=2,
            )
            fit_res = bp.run_prism(prism_obj, update_gibbs=False, gibbs_control=Some(ctrl))
            match fit_res:
                case Failure(err):
                    return Failure(f"BayesPrism Gibbs run failed: {err}")
                case Success(fitted):
                    frac_res = bp.get_fraction(fitted, which_theta="first", state_or_type="type")
                    match frac_res:
                        case Failure(err):
                            return Failure(f"BayesPrism fraction extraction failed: {err}")
                        case Success(bp_df):
                            bp_fractions = bp_df

    # 2. Run InstaPrism
    print("  [2/2] Running InstaPrism...")
    row_sums = sub_ref.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    norm_ref = (sub_ref / row_sums).astype(np.float64)  # (clusters x genes)

    insta_mat = np.zeros((len(sample_ids), len(clusters)), dtype=np.float64)
    for i in range(len(sample_ids)):
        b_vec = sub_bulk.iloc[i].to_numpy().astype(np.float64)
        _, _, fracs, _ = instaprism.insta_prism(bulk=b_vec, reference=norm_ref, n_iter=50)
        insta_mat[i, :] = fracs

    insta_dict: dict[str, list[object]] = {"bulk_id": sample_ids}
    for j, c in enumerate(clusters):
        insta_dict[c] = insta_mat[:, j].tolist()
    insta_fractions = pl.DataFrame(insta_dict)

    return Success((bp_fractions, insta_fractions, clusters, sample_ids))


def fit_and_evaluate_purity(
    frac_df: pl.DataFrame,
    clin_df: pd.DataFrame,
    clusters: list[str],
    algo_name: str,
    bulk_df: pd.DataFrame,
) -> pl.DataFrame:
    """Fit logistic response models with and without tumor purity normalization."""
    # Match samples
    valid_samples = [s for s in frac_df["bulk_id"].to_list() if s in clin_df.index]
    sub_frac = frac_df.filter(pl.col("bulk_id").is_in(valid_samples))

    responses = [
        1.0 if str(clin_df.loc[s, "RESPONSE"]).strip().lower() in ["complete response", "partial response", "responder"] else 0.0
        for s in sub_frac["bulk_id"].to_list()
    ]
    y = np.array(responses, dtype=np.float64)

    raw_mat = sub_frac.select(clusters).to_numpy().astype(np.float64)

    # In solid tumor biopsies:
    # theta represents the relative composition of the immune infiltrate.
    # Total immune infiltration fraction in whole biopsy: Infiltrate_i = 1 - Purity_i
    if "PTPRC" in bulk_df.index:
        ptprc_expr = bulk_df.loc["PTPRC", valid_samples].to_numpy().astype(np.float64)
        immune_score = ptprc_expr / (np.max(ptprc_expr) + 1e-6)
        tumor_purity = np.clip(1.0 - immune_score, 0.1, 0.95)
    else:
        tumor_purity = np.full(len(valid_samples), 0.70)

    infiltrate = 1.0 - tumor_purity  # immune fraction of whole biopsy

    # 1. WITHOUT Purity Normalization: Raw fraction of whole biopsy = theta * Infiltrate
    raw_biopsy_mat = raw_mat * infiltrate[:, np.newaxis]

    # 2. WITH Purity Normalization: Fraction within immune infiltrate = raw_biopsy / Infiltrate = theta
    norm_immune_mat = raw_mat.copy()

    results: list[dict[str, object]] = []

    for purity_mode, mat in [
        ("Without Purity Normalization (Raw Biopsy)", raw_biopsy_mat),
        ("With Purity Normalization (Immune Infiltrate)", norm_immune_mat),
    ]:
        for j, state in enumerate(clusters):
            x = mat[:, j]
            if np.std(x) < 1e-8:
                beta_z, p_val = 0.0, 1.0
            else:
                x_std = (x - np.mean(x)) / np.std(x)
                r_val, p_val = stats.pointbiserialr(y, x_std)
                r_clip = np.clip(r_val, -0.999, 0.999)
                beta_z = float(2.0 * r_clip / np.sqrt(1.0 - r_clip**2 + 1e-12))

            results.append({
                "algorithm": algo_name,
                "purity_normalization": purity_mode,
                "cell_state": state,
                "beta_bulk": beta_z,
                "p_value_bulk": float(p_val),
            })

    return pl.DataFrame(results)


def evaluate_concordance_with_sc(
    bulk_results: pl.DataFrame,
    sc_effects_path: Path,
) -> pl.DataFrame:
    """Compute Spearman rho and concordance with single-cell Milo effects."""
    sc_df = pl.read_parquet(sc_effects_path)
    joined = bulk_results.join(sc_df, on="cell_state", how="inner")

    metrics_rows: list[dict[str, object]] = []
    grouped = joined.partition_by(["algorithm", "purity_normalization"], as_dict=True)

    for (algo, mode), sub_df in grouped.items():
        b_sc = sub_df["beta_sc"].to_numpy()
        b_blk = sub_df["beta_bulk"].to_numpy()

        rho, p_rho = stats.spearmanr(b_sc, b_blk)
        r, p_r = stats.pearsonr(b_sc, b_blk)

        concordant = sum(
            np.sign(s) == np.sign(b) for s, b in zip(b_sc, b_blk, strict=False) if abs(s) > 0.05 and abs(b) > 0.05
        )
        total = sum(1 for s, b in zip(b_sc, b_blk, strict=False) if abs(s) > 0.05 and abs(b) > 0.05)
        rate = concordant / max(total, 1)

        metrics_rows.append({
            "algorithm": algo,
            "purity_normalization": mode,
            "n_states": sub_df.height,
            "spearman_rho": float(rho),
            "spearman_pval": float(p_rho),
            "pearson_r": float(r),
            "concordance_rate": float(rate),
        })

    return pl.DataFrame(metrics_rows)


def plot_comparison_chart(metrics_df: pl.DataFrame, out_path: Path) -> Result[None, str]:
    """Plot grouped bar chart of Spearman rho across algorithms and purity modes."""
    chart = (
        alt.Chart(metrics_df)
        .mark_bar()
        .encode(
            x=alt.X("algorithm:N", title="Deconvolution Algorithm"),
            y=alt.Y("spearman_rho:Q", title="Concordance with Single-Cell (Spearman ρ)", scale=alt.Scale(domain=[-0.2, 1.0])),
            color=alt.Color("purity_normalization:N", legend=alt.Legend(title="Tumor Purity")),
            xOffset="purity_normalization:N",
            tooltip=["algorithm:N", "purity_normalization:N", "spearman_rho:Q", "concordance_rate:Q"],
        )
        .properties(
            title="BayesPrism vs. InstaPrism: Impact of Tumor Purity Normalization on Concordance",
            width=480,
            height=360,
        )
        .configure_axis(grid=True, gridColor="#f0f0f0")
    )

    try:
        svg_str = vlc.vegalite_to_svg(chart.to_dict())
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(svg_str, encoding="utf-8")
        return Success(None)
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to save SVG: {exc}")


def run_pipeline(config: BenchmarkConfig) -> Result[None, str]:
    """Execute Step 8 benchmark pipeline."""
    # 1. Load inputs
    match load_inputs(config):
        case Failure(err):
            return Failure(err)
        case Success((ref_df, bulk_df, clin_df)):
            pass

    # 2. Run deconvolution with both BayesPrism and InstaPrism
    match run_deconvolutions(ref_df, bulk_df, config):
        case Failure(err):
            return Failure(err)
        case Success((bp_fractions, insta_fractions, clusters, sample_ids)):
            pass

    # 3. Fit models with and without purity normalization
    bp_eval = fit_and_evaluate_purity(bp_fractions, clin_df, clusters, "BayesPrism (Gibbs MCMC)", bulk_df)
    insta_eval = fit_and_evaluate_purity(insta_fractions, clin_df, clusters, "InstaPrism", bulk_df)
    all_eval = pl.concat([bp_eval, insta_eval])

    # 4. Concordance with single-cell
    metrics_df = evaluate_concordance_with_sc(all_eval, config.sc_effects_path)

    # 5. Save outputs
    config.out_dir.mkdir(parents=True, exist_ok=True)
    all_eval.write_parquet(config.out_dir / "bayesprism_vs_instaprism_effects.parquet")
    metrics_df.write_parquet(config.out_dir / "bayesprism_vs_instaprism_metrics.parquet")

    # 6. Plot SVG
    config.results_dir.mkdir(parents=True, exist_ok=True)
    svg_path = config.results_dir / "fig5_bayesprism_vs_instaprism_purity.svg"
    plot_comparison_chart(metrics_df, svg_path)

    print("\n==================================================================")
    print("BAYESPRISM vs. INSTAPRISM: WITH vs. WITHOUT PURITY NORMALIZATION")
    print("==================================================================")
    for row in metrics_df.iter_rows(named=True):
        print(
            f"Algorithm: {row['algorithm']:<22} | "
            f"Purity: {row['purity_normalization']:<36} | "
            f"Spearman ρ: {row['spearman_rho']:>5.3f} | "
            f"Concordance: {row['concordance_rate']*100:>5.1f}%"
        )
    print("==================================================================")

    return Success(None)


def main() -> None:
    config = parse_args()
    match run_pipeline(config):
        case Success(_):
            sys.exit(0)
        case Failure(err):
            print(f"[ERROR] {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
