"""Unit and property-based tests for Bagaev MFP 29-signature ssGSEA and median/MAD scaling."""

from __future__ import annotations

from pathlib import Path
import anndata as ad  # type: ignore[import-untyped]
from hypothesis import given, settings, strategies as st
import numpy as np
import pandas as pd
import polars as pl
from returns.maybe import Nothing, Some
from returns.result import Failure, Success
import scipy.sparse as sp  # type: ignore[import-untyped]

from tme_datasets import (
    BAGAEV_MFP_29_SIGNATURES,
    BAGAEV_MFP_FUNCTIONAL_GROUPS,
    BAGAEV_MFP_SIGNATURE_ORDER,
    GeneSet,
    GeneSetCollection,
    TumorPseudobulkSpec,
    compute_ssgsea_matrix,
    deduplicate_patient_profiles,
    get_bagaev_mfp_collection,
    load_ensembl_symbol_map,
    load_geneset_collection,
    load_panmelanoma_bagaev_labels,
    load_tcga_bagaev_reference,
    match_iatlas_to_panmelanoma_labels,
    median_mad_scale_signatures,
    score_geneset_ssgsea,
    stream_tumor_pseudobulk_h5ad,
)


def test_bagaev_mfp_29_collection_completeness_and_order() -> None:
    """Verify all 29 Bagaev Fges signatures are present, ordered, and partitioned into 4 groups."""
    col = get_bagaev_mfp_collection()
    assert col.id == "bagaev_mfp_29"
    assert len(BAGAEV_MFP_SIGNATURE_ORDER) == 29
    assert len(BAGAEV_MFP_29_SIGNATURES) == 29
    assert tuple(col.gene_sets.keys()) == BAGAEV_MFP_SIGNATURE_ORDER

    grouped_sigs = [s for grp in BAGAEV_MFP_FUNCTIONAL_GROUPS.values() for s in grp]
    assert len(grouped_sigs) == 29
    assert set(grouped_sigs) == set(BAGAEV_MFP_SIGNATURE_ORDER)

    loaded = load_geneset_collection("bagaev_mfp_29")
    assert isinstance(loaded, Success)
    assert tuple(loaded.unwrap().gene_sets.keys()) == BAGAEV_MFP_SIGNATURE_ORDER


def test_ssgsea_matches_reference_formula_exactly() -> None:
    """Verify compute_ssgsea_matrix reproduces the exact Bagaev rank-exponent formula."""
    rng = np.random.default_rng(123)
    genes = ("HLA-A", "HLA-B", "CD8A", "GZMB", "COL1A1", "FAP", "MKI67", "VEGFA")
    mat = rng.uniform(0.0, 15.0, size=(5, len(genes)))
    samples = tuple(f"S{i}" for i in range(5))
    collection = GeneSetCollection(
        id="mini",
        name="Mini",
        gene_sets={
            "Eff": GeneSet(id="Eff", name="Eff", genes=("CD8A", "GZMB")),
            "CAF": GeneSet(id="CAF", name="CAF", genes=("COL1A1", "FAP")),
        },
    )

    res = compute_ssgsea_matrix(mat, samples, genes, collection)
    assert isinstance(res, Success)
    df = res.unwrap()
    assert df.shape == (5, 3)

    # Compare against explicit pandas rank formula from Bagaev/BostonGene MFP
    df_pd = pd.DataFrame(mat, index=list(samples), columns=list(genes))
    ranks_pd = df_pd.T.rank(method="max", na_option="bottom")
    for sig_name, sig_genes in (("Eff", ["CD8A", "GZMB"]), ("CAF", ["COL1A1", "FAP"])):
        sr = ranks_pd.loc[sig_genes]
        expected = (sr**1.25).sum(axis=0) / (sr**0.25).sum(axis=0) - (len(genes) - len(sig_genes) + 1) / 2.0
        np.testing.assert_allclose(df[sig_name].to_numpy(), expected.to_numpy(), rtol=1e-10, atol=1e-10)


@given(
    scale=st.floats(min_value=0.1, max_value=10.0, allow_nan=False, allow_infinity=False),
    shift=st.floats(min_value=0.0, max_value=50.0, allow_nan=False, allow_infinity=False),
)
@settings(max_examples=25)
def test_ssgsea_strictly_invariant_under_monotonic_transforms(scale: float, shift: float) -> None:
    """Property: positive affine & log1p monotonic transforms within samples preserve ssGSEA identically."""
    base = np.array(
        [
            [0.5, 2.1, 8.4, 15.2, 0.0, 4.3],
            [12.0, 1.0, 0.2, 5.5, 9.9, 3.1],
        ],
        dtype=np.float64,
    )
    genes = ("CD8A", "GZMB", "PRF1", "COL1A1", "FAP", "VEGFA")
    samples = ("P1", "P2")
    col = GeneSetCollection(
        id="prop",
        name="Prop",
        gene_sets={"T": GeneSet(id="T", name="T", genes=("CD8A", "GZMB", "PRF1"))},
    )

    r1 = compute_ssgsea_matrix(base, samples, genes, col)
    r2 = compute_ssgsea_matrix(np.log1p(base * scale + shift), samples, genes, col)
    assert isinstance(r1, Success) and isinstance(r2, Success)
    np.testing.assert_allclose(r1.unwrap()["T"].to_numpy(), r2.unwrap()["T"].to_numpy(), atol=1e-12)


def test_score_geneset_ssgsea_resolves_ensembl_var_gene_name() -> None:
    """Verify score_geneset_ssgsea maps Ensembl var_names via adata.var['gene_name']."""
    mat = np.array([[10.0, 9.0, 1.0, 0.5], [0.5, 1.0, 11.0, 12.0]], dtype=np.float32)
    adata = ad.AnnData(
        X=mat,
        obs=pd.DataFrame(index=["S1", "S2"]),
        var=pd.DataFrame(
            {"gene_name": ["CD8A", "GZMB", "COL1A1", "FAP"]},
            index=["ENSG00000153563", "ENSG00000100453", "ENSG00000108821", "ENSG00000078098"],
        ),
    )
    col = GeneSetCollection(
        id="t",
        name="t",
        gene_sets={
            "Eff": GeneSet(id="Eff", name="Eff", genes=("CD8A", "GZMB")),
            "CAF": GeneSet(id="CAF", name="CAF", genes=("COL1A1", "FAP")),
        },
    )
    res = score_geneset_ssgsea(adata, col)
    assert isinstance(res, Success)
    df = res.unwrap()
    assert df["Eff"][0] > df["Eff"][1]
    assert df["CAF"][1] > df["CAF"][0]


def test_median_mad_scale_signatures_global_and_grouped() -> None:
    """Verify global and cohort-grouped median-centering + MAD scaling and clipping."""
    df = pl.DataFrame(
        {
            "sample_id": ["A1", "A2", "A3", "B1", "B2", "B3"],
            "cohort": ["C1", "C1", "C1", "C2", "C2", "C2"],
            "CAF": [1.0, 2.0, 3.0, 101.0, 102.0, 103.0],
        }
    )
    res_grp = median_mad_scale_signatures(df, ("CAF",), group_col=Some("cohort"), clip=Some(4.0))
    assert isinstance(res_grp, Success)
    scaled_grp = res_grp.unwrap()["CAF"].to_numpy()
    # Both cohorts have identical relative spacing [-1.5, 0.0, +1.5] after cohort-wise median/MAD scaling
    np.testing.assert_allclose(scaled_grp[:3], scaled_grp[3:], atol=1e-10)
    assert np.isclose(scaled_grp[1], 0.0)

    err_res = median_mad_scale_signatures(df, ("MISSING_COL",), group_col=Nothing)
    assert isinstance(err_res, Failure)


def test_load_tcga_bagaev_reference_and_panmelanoma_timepoint_preservation(tmp_path: Path) -> None:
    """Verify TCGA reference loader and non-colliding Pre/On Pan-Melanoma matching."""
    tcga_dir = tmp_path / "Bagaev_TCGA"
    tcga_dir.mkdir(parents=True)
    (tcga_dir / "annotations.tsv").write_text(
        "\tTCGA_project\tHISTOLOGICAL_SUBTYPE\tMFP\n"
        "TCGA-01\tSKCM\tCutaneous\tIE\n"
        "TCGA-02\tSKCM\tCutaneous\tF\n",
        encoding="utf-8",
    )
    (tcga_dir / "signature.tsv").write_text(
        "\tTCGA-01\tTCGA-02\nEffector_cells\t1.5\t-1.2\nCAF\t-0.9\t2.1\n",
        encoding="utf-8",
    )
    ref_res = load_tcga_bagaev_reference(Some(tmp_path))
    assert isinstance(ref_res, Success)
    ref_df = ref_res.unwrap()
    assert ref_df.shape == (2, 6)
    assert ref_df.filter(pl.col("sample_id") == "TCGA-01")["MFP"][0] == "IE"

    pan_dir = tmp_path / "Bagaev_PanMelanoma"
    pan_dir.mkdir(parents=True)
    # Same Patient (Riaz_Pt1) switches MFP across Pre vs On — must NOT overwrite Pre with On
    (pan_dir / "annotations.tsv").write_text(
        "Sample\tCohort\tMFP\tPatient\tResponse\tTime\n"
        "SRR1\tRiaz\tIE\tRiaz_Pt1\tR\tPre\n"
        "SRR2\tRiaz\tD\tRiaz_Pt1\tR\tOn\n",
        encoding="utf-8",
    )
    pan_res = load_panmelanoma_bagaev_labels(Some(tmp_path))
    assert isinstance(pan_res, Success)
    iatlas_df = pl.DataFrame(
        {
            "cohort": ["Riaz-iAtlas", "Riaz-iAtlas"],
            "raw_sample_id": ["Pt1_PRE", "Pt1_ON"],
            "patient_id": ["Pt1", "Pt1"],
            "sample_treatment": ["Pre", "On"],
        }
    )
    matched = match_iatlas_to_panmelanoma_labels(iatlas_df, pan_res.unwrap())
    assert matched["published_bagaev_mfp"].to_list() == ["IE", "D"]


def test_stream_tumor_pseudobulk_excludes_blood_resolves_ensg_and_inverts_log_scales(tmp_path: Path) -> None:
    """Verify tumor-only gating, Ensembl parquet resolution, ln vs log2 base inversion & dedup."""
    # Write synthetic release 111 Ensembl map parquet
    ens_dir = tmp_path / "ensembl"
    ens_dir.mkdir()
    pl.DataFrame(
        {
            "gene_id": ["ENSG00000153563", "ENSG00000108821"],
            "gene_name": ["CD8A", "COL1A1"],
        }
    ).write_parquet(ens_dir / "gene_mapping_cache_release_111.parquet")
    ens_map = load_ensembl_symbol_map(ens_dir)

    # Cohort 1: CELLxGENE_05a8c945 style with blood + normal + tumor + metastasis and ln(CP10k + 1) X
    # Donor D1 has a blood sample with high CD8A and a tumor sample with high COL1A1; blood must be excluded!
    n_blood = 30
    n_tum = 30
    x_ln = np.vstack(
        [
            np.tile(np.array([[np.log1p(9000.0), np.log1p(1000.0)]]), (n_blood, 1)),
            np.tile(np.array([[np.log1p(1000.0), np.log1p(9000.0)]]), (n_tum, 1)),
        ]
    ).astype(np.float32)
    ad_crc = ad.AnnData(
        X=sp.csr_matrix(x_ln),
        obs=pd.DataFrame(
            {
                "sample_id": ["D1_Blood"] * n_blood + ["D1_Tumor"] * n_tum,
                "donor_id": ["D1"] * (n_blood + n_tum),
                "sample_type": ["blood"] * n_blood + ["tumor"] * n_tum,
                "treatment_status": ["naive"] * (n_blood + n_tum),
                "clinical_response": ["responder"] * (n_blood + n_tum),
            },
            index=[f"c_{i}" for i in range(n_blood + n_tum)],
        ),
        var=pd.DataFrame(index=["ENSG00000153563.10", "ENSG00000108821"]),
    )
    crc_path = tmp_path / "crc.h5ad"
    ad_crc.write_h5ad(crc_path)

    res_crc = stream_tumor_pseudobulk_h5ad(
        crc_path,
        TumorPseudobulkSpec(cohort_id="CELLxGENE_05a8c945", sample_col="sample_id", min_cells_per_sample=20),
        ens_map,
    )
    assert isinstance(res_crc, Success)
    out_crc = res_crc.unwrap()
    assert out_crc.scale_mode == "ln_cp10k"
    assert out_crc.sample_ids == ("D1_Tumor",)
    assert out_crc.gene_symbols == ("CD8A", "COL1A1")
    # Reconstructed log2(CPM+1) must preserve COL1A1:CD8A = 9:1 -> COL1A1 ~ 900k CPM, CD8A ~ 100k CPM
    assert np.isclose(out_crc.log2_cpm[0, 1], np.log2(900000.0 + 1.0), atol=1e-2)
    assert np.isclose(out_crc.log2_cpm[0, 0], np.log2(100000.0 + 1.0), atol=1e-2)

    # Verify log2(TPM + 1) Smart-seq2 detection (max_nz > 12.5) and GSE120575 patient deduplication
    x_log2 = np.vstack(
        [
            np.tile(np.array([[np.log2(65536.0 + 1.0), np.log2(1024.0 + 1.0)]]), (25, 1)),
            np.tile(np.array([[np.log2(1024.0 + 1.0), np.log2(65536.0 + 1.0)]]), (35, 1)),
        ]
    ).astype(np.float32)
    ad_ss2 = ad.AnnData(
        X=sp.csr_matrix(x_log2),
        obs=pd.DataFrame(
            {
                "patient_geo": ["Pre_P1"] * 25 + ["Post_P1"] * 35,
                "response": ["Responder"] * 60,
            },
            index=[f"s_{i}" for i in range(60)],
        ),
        var=pd.DataFrame(index=["CD8A", "COL1A1"]),
    )
    ss2_path = tmp_path / "gse120575.h5ad"
    ad_ss2.write_h5ad(ss2_path)
    res_ss2 = stream_tumor_pseudobulk_h5ad(
        ss2_path,
        TumorPseudobulkSpec(cohort_id="GSE120575", sample_col="patient_geo", min_cells_per_sample=20),
        ens_map,
    )
    assert isinstance(res_ss2, Success)
    out_ss2 = res_ss2.unwrap()
    assert out_ss2.scale_mode == "log2_tpm"
    dedup_df = deduplicate_patient_profiles(
        out_ss2.meta_df.with_columns(pl.lit("GSE120575").alias("cohort")),
        mode="strict_pre",
        weight_col=Some("n_cells"),
    )
    assert dedup_df.height == 1
    assert dedup_df["sample_id"][0] == "Pre_P1"
    assert dedup_df["patient_id"][0] == "P1"
