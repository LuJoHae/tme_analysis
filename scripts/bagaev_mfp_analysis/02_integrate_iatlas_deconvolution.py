#!/usr/bin/env python3
"""Step 02: Integrate Bagaev TME Subtypes with Consensus-Meta42 & SC9 Deconvolution Predictions.

Analyzes:
1. Enrichment of deconvolved Consensus-Meta42 cell states (M01..M42, annotated with lineage &
   top marker genes) and SC9 multi-resolution states across Bagaev IE, IE/F, F, and D subtypes.
2. Alignment between Bagaev's 4 subtypes (especially IE vs IE/F fibrotic suppression) and our
   Consensus-Meta42 2x2 Stromal-Gating / Epistasis regimes (`oof_9a_9b`).
3. Deconvolution response prediction distributions across Bagaev subtypes AND within-subtype
   responder discrimination (ROC-AUC / PR-AUC inside IE, IE/F, F, and D individually).
4. Cross-validated head-to-head comparison: Bagaev 4-Subtype categorical vs Bagaev 29-Fges
   continuous vs Consensus-Meta42 Deconvolution OOF vs Combined (Bagaev + Deconvolution).
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import polars as pl
from returns.result import Failure, Success
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.preprocessing import StandardScaler

from tme_datasets import BAGAEV_MFP_SIGNATURE_ORDER, get_data_paths
from tme_response import (
    CANONICAL_SUBTYPE_ORDER,
    analyze_subtype_deconvolution_enrichment,
    evaluate_within_subtype_deconv_predictions,
)


def _evaluate_cv_models(joined_df: pl.DataFrame) -> pl.DataFrame:
    """Compare 5x5 repeated stratified CV ROC-AUC & PR-AUC across Bagaev, Deconv, and Combined models."""
    clean = joined_df.filter(
        pl.col("response_binary").cast(pl.Float64).is_in([0.0, 1.0])
        & pl.col("oof_deconv_net_z").is_not_null()
        & pl.col("bagaev_subtype").is_not_null()
    )
    y = clean["response_binary"].to_numpy().astype(np.int64)

    # Feature matrices
    subtypes = clean["bagaev_subtype"].cast(pl.String).to_list()
    x_onehot = np.column_stack(
        [[1.0 if s == st_name else 0.0 for s in subtypes] for st_name in CANONICAL_SUBTYPE_ORDER]
    )
    fges_cols = [c for c in BAGAEV_MFP_SIGNATURE_ORDER if c in clean.columns]
    x_fges = clean.select(fges_cols).to_numpy().astype(np.float64)
    meta_cols = [f"M{i:02d}" for i in range(1, 43) if f"M{i:02d}" in clean.columns]
    x_meta = np.arcsin(np.sqrt(np.clip(clean.select(meta_cols).to_numpy().astype(np.float64), 0.0, 1.0)))
    oof_pred = clean["oof_deconv_net_z"].to_numpy().astype(np.float64)

    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=42)
    model_specs = (
        ("Bagaev_4Subtype_Categorical", x_onehot),
        ("Bagaev_29Fges_L2_Logistic", x_fges),
        ("Deconv_Meta42_Fractions_L2", x_meta),
        ("Combined_Bagaev29_Plus_Meta42", np.hstack([x_fges, x_meta])),
    )

    rows: list[dict[str, str | float]] = [
        {
            "model": "Published_OOF_Meta42_StromalGating_3c",
            "roc_auc_mean": float(roc_auc_score(y, oof_pred)),
            "roc_auc_sd": 0.0,
            "pr_auc_mean": float(average_precision_score(y, oof_pred)),
            "pr_auc_sd": 0.0,
        }
    ]

    for name, x_mat in model_specs:
        aucs: list[float] = []
        praucs: list[float] = []
        for tr_idx, te_idx in rskf.split(x_mat, y):
            scaler = StandardScaler()
            x_tr = scaler.fit_transform(x_mat[tr_idx])
            x_te = scaler.transform(x_mat[te_idx])
            clf = LogisticRegression(C=0.25, max_iter=1000, class_weight="balanced")
            clf.fit(x_tr, y[tr_idx])
            prob = clf.predict_proba(x_te)[:, 1]
            aucs.append(float(roc_auc_score(y[te_idx], prob)))
            praucs.append(float(average_precision_score(y[te_idx], prob)))
        rows.append(
            {
                "model": name,
                "roc_auc_mean": float(np.mean(aucs)),
                "roc_auc_sd": float(np.std(aucs)),
                "pr_auc_mean": float(np.mean(praucs)),
                "pr_auc_sd": float(np.std(praucs)),
            }
        )

    # Also test per-cohort mean CV ROC-AUC for cohorts with >=25 samples and >=5 responders
    cohort_rows: list[dict[str, str | float]] = []
    for grp in clean.partition_by("cohort", maintain_order=True):
        cid = str(grp["cohort"][0])
        yc = grp["response_binary"].to_numpy().astype(np.int64)
        if grp.height < 25 or int(np.sum(yc == 1)) < 5 or int(np.sum(yc == 0)) < 5:
            continue
        p_oof = grp["oof_deconv_net_z"].to_numpy().astype(np.float64)
        p_ie_prob = (grp["prob_IE"] + 0.5 * grp["prob_IEF"]).to_numpy().astype(np.float64)
        cohort_rows.append(
            {
                "cohort": cid,
                "n_samples": float(grp.height),
                "auc_bagaev_ie_score": float(roc_auc_score(yc, p_ie_prob)),
                "auc_deconv_oof": float(roc_auc_score(yc, p_oof)),
                "auc_combined_rank": float(
                    roc_auc_score(
                        yc,
                        0.5 * (p_ie_prob - np.mean(p_ie_prob)) / max(1e-8, np.std(p_ie_prob))
                        + 0.5 * (p_oof - np.mean(p_oof)) / max(1e-8, np.std(p_oof)),
                    )
                ),
            }
        )

    out_df = pl.DataFrame(rows)
    if cohort_rows:
        out_df = out_df.with_columns(pl.lit("PanCancer_CV").alias("scope"))
    return out_df


def run_step02(
    bagaev_dir: Path,
    bench_3c_dir: Path,
    sc9_deconv_dir: Path,
) -> None:
    """Execute Step 02 integration of iAtlas Bagaev subtypes with deconvolution outputs."""
    subtypes_path = bagaev_dir / "iatlas_bagaev_sample_subtypes.parquet"
    if not subtypes_path.is_file():
        raise FileNotFoundError(f"Step 01 output not found at {subtypes_path}")
    subtypes_df = pl.read_parquet(subtypes_path)

    frac_path = bench_3c_dir / "fractions" / "fractions_Consensus-Meta42_on_PanCancer-Combined.parquet"
    oof_path = bench_3c_dir / "oof_9a_9b_patient_landscape_and_curves.parquet"
    hub_path = bench_3c_dir / "consensus_metastate_all42_hub_summary.parquet"
    deg_path = bench_3c_dir / "consensus_metastate_all42_top_degs.parquet"

    if not frac_path.is_file() or not oof_path.is_file():
        raise FileNotFoundError(f"Missing 3c deconvolution outputs in {bench_3c_dir}")

    frac_df = pl.read_parquet(frac_path)
    meta_cols = [f"M{i:02d}" for i in range(1, 43) if f"M{i:02d}" in frac_df.columns]
    oof_patients = (
        pl.read_parquet(oof_path)
        .filter(pl.col("record_type") == "patient")
        .with_columns((pl.col("x_z") - pl.col("y_z")).alias("oof_deconv_net_z"))
        .select(["sample_id", "regime_group", "stratum_label", "x_z", "y_z", "oof_deconv_net_z"])
    )

    joined = (
        subtypes_df.join(frac_df.select(["sample_id", *meta_cols]), on="sample_id", how="inner")
        .join(oof_patients, on="sample_id", how="left")
    )
    joined.write_parquet(bagaev_dir / "iatlas_bagaev_joined_meta42_oof.parquet")
    print(f"Joined {joined.height} evaluable iAtlas patients with Consensus-Meta42 fractions & OOF predictions.")

    # 1. Meta-state enrichment across IE, IE/F, F, D
    enr_res = analyze_subtype_deconvolution_enrichment(joined, meta_cols, subtype_col="bagaev_subtype")
    if isinstance(enr_res, Failure):
        raise RuntimeError(f"Meta42 enrichment failed: {enr_res.failure()}")
    enr_df = enr_res.unwrap()

    if hub_path.is_file():
        hub_df = pl.read_parquet(hub_path)
        id_col = next((c for c in ("meta_state", "state_id", "metastate_id") if c in hub_df.columns), None)
        if id_col is not None:
            keep_cols = [c for c in hub_df.columns if c != id_col and c not in enr_df.columns][:6]
            enr_df = enr_df.join(
                hub_df.select([pl.col(id_col).alias("cell_state"), *keep_cols]),
                on="cell_state",
                how="left",
            )
    if deg_path.is_file():
        deg_df = pl.read_parquet(deg_path)
        st_c = next((c for c in ("meta_state", "state_id", "cell_state") if c in deg_df.columns), None)
        gn_c = next((c for c in ("gene_symbol", "gene_name", "gene") if c in deg_df.columns), None)
        if st_c is not None and gn_c is not None:
            top_markers = (
                deg_df.group_by(st_c)
                .agg(pl.col(gn_c).head(5).str.join(", ").alias("top_deg_markers"))
                .rename({st_c: "cell_state"})
            )
            enr_df = enr_df.join(top_markers, on="cell_state", how="left")

    enr_df.write_parquet(bagaev_dir / "iatlas_meta42_enrichment_by_bagaev_subtype.parquet")
    enr_df.write_csv(bagaev_dir / "iatlas_meta42_enrichment_by_bagaev_subtype.tsv", separator="\t")

    # 2. Cross-tabulation of Bagaev subtypes vs 3c Stromal-Gating / Epistasis Regimes
    if "regime_group" in joined.columns:
        regime_cross = (
            joined.filter(pl.col("regime_group").is_not_null())
            .group_by(["bagaev_subtype", "regime_group"])
            .agg(
                [
                    pl.len().alias("n_patients"),
                    pl.col("response_binary").mean().alias("observed_orr"),
                    pl.col("x_z").mean().alias("mean_oof_responder_z"),
                    pl.col("y_z").mean().alias("mean_oof_exclusion_z"),
                    pl.col("oof_deconv_net_z").mean().alias("mean_oof_deconv_net_z"),
                ]
            )
            .sort(["bagaev_subtype", "regime_group"])
        )
        regime_cross.write_parquet(bagaev_dir / "iatlas_bagaev_vs_stromal_gating_regimes.parquet")
        regime_cross.write_csv(bagaev_dir / "iatlas_bagaev_vs_stromal_gating_regimes.tsv", separator="\t")

    # 3. Within-subtype OOF prediction discrimination (Does deconvolution stratify inside IE and IE/F?)
    within_res = evaluate_within_subtype_deconv_predictions(
        joined, score_col="oof_deconv_net_z", subtype_col="bagaev_subtype", response_col="response_binary"
    )
    if isinstance(within_res, Failure):
        raise RuntimeError(f"Within-subtype evaluation failed: {within_res.failure()}")
    within_df = within_res.unwrap()
    within_df.write_parquet(bagaev_dir / "iatlas_within_subtype_deconv_auc.parquet")
    within_df.write_csv(bagaev_dir / "iatlas_within_subtype_deconv_auc.tsv", separator="\t")

    # 4. Cross-validated model benchmark (Bagaev vs Meta42 Deconv vs Combined)
    cv_bench = _evaluate_cv_models(joined)
    cv_bench.write_parquet(bagaev_dir / "iatlas_bagaev_vs_deconv_cv_benchmark.parquet")
    cv_bench.write_csv(bagaev_dir / "iatlas_bagaev_vs_deconv_cv_benchmark.tsv", separator="\t")

    # 5. Also check SC9 multi-resolution combined reference fractions (res0.5) if present
    sc9_frac_file = sc9_deconv_dir / "deconv_fractions_all_9_combined_res0.5.parquet"
    if sc9_frac_file.is_file():
        sc9_df = pl.read_parquet(sc9_frac_file)
        if "sample_id" in sc9_df.columns:
            sc9_state_cols = [
                c for c in sc9_df.columns
                if c not in ("sample_id", "cohort", "cohort_id", "response", "response_binary", "dataset_id")
                and sc9_df[c].dtype in (pl.Float32, pl.Float64)
            ]
            # Match either prefixed or unprefixed sample_ids
            sc9_join = subtypes_df.join(sc9_df.select(["sample_id", *sc9_state_cols]), on="sample_id", how="inner")
            if sc9_join.height < 50 and "raw_sample_id" in subtypes_df.columns:
                sc9_join = subtypes_df.join(
                    sc9_df.select([pl.col("sample_id").alias("raw_sample_id"), *sc9_state_cols]),
                    on="raw_sample_id",
                    how="inner",
                )
            if sc9_join.height >= 20 and sc9_state_cols:
                sc9_enr = analyze_subtype_deconvolution_enrichment(sc9_join, sc9_state_cols)
                if isinstance(sc9_enr, Success):
                    sc9_enr.unwrap().write_parquet(bagaev_dir / "iatlas_sc9_res05_enrichment_by_bagaev_subtype.parquet")

    print(f"Step 02 complete. Saved deconvolution integration tables to {bagaev_dir}")


def main() -> None:
    paths = get_data_paths()
    local_3c = Path("output/manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark")
    default_3c = local_3c if (local_3c / "oof_9a_9b_patient_landscape_and_curves.parquet").is_file() else (paths.output_dir / "manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark")
    local_sc9 = Path("output/sc9_milo_deconv_validation/02_deconvolution")
    default_sc9 = local_sc9 if local_sc9.is_dir() else (paths.output_dir / "sc9_milo_deconv_validation/02_deconvolution")

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bagaev-dir", type=Path, default=paths.output_dir / "bagaev_mfp_analysis")
    parser.add_argument("--bench-3c-dir", type=Path, default=default_3c)
    parser.add_argument("--sc9-deconv-dir", type=Path, default=default_sc9)
    args = parser.parse_args()
    run_step02(bagaev_dir=args.bagaev_dir, bench_3c_dir=args.bench_3c_dir, sc9_deconv_dir=args.sc9_deconv_dir)


if __name__ == "__main__":
    main()
