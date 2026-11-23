"""Unit and property-based tests for Bagaev MFP 29-signature ssGSEA and median/MAD scaling."""

from __future__ import annotations

import anndata as ad  # type: ignore[import-untyped]
from hypothesis import given, settings, strategies as st
import numpy as np
import pandas as pd
import polars as pl
from returns.maybe import Nothing, Some
from returns.result import Failure, Success

from tme_datasets import (
    BAGAEV_MFP_29_SIGNATURES,
    BAGAEV_MFP_FUNCTIONAL_GROUPS,
    BAGAEV_MFP_SIGNATURE_ORDER,
    GeneSet,
    GeneSetCollection,
    compute_ssgsea_matrix,
    get_bagaev_mfp_collection,
    load_geneset_collection,
    median_mad_scale_signatures,
    score_geneset_ssgsea,
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
