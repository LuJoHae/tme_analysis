"""Unit tests for gene sets, collections, and signature scoring."""

from pathlib import Path
import anndata as ad
from returns.result import Success
from tme_datasets.genesets.collections import get_tme_major_lineage_collection
from tme_datasets.genesets.overlap import compute_geneset_overlap
from tme_datasets.genesets.parser import export_gmt, parse_gmt
from tme_datasets.genesets.scoring import score_geneset_auc, score_geneset_zscore


def test_tme_major_lineage_collection() -> None:
    coll = get_tme_major_lineage_collection()
    assert "T_NK" in coll.gene_sets
    assert "CD8A" in coll.gene_sets["T_NK"].genes


def test_gmt_export_and_parse(tmp_path: Path) -> None:
    coll = get_tme_major_lineage_collection()
    gmt_file = tmp_path / "tme_markers.gmt"

    exp_res = export_gmt(coll, gmt_file)
    assert isinstance(exp_res, Success)
    assert gmt_file.exists()

    parse_res = parse_gmt(gmt_file, collection_id="parsed_tme")
    assert isinstance(parse_res, Success)
    parsed = parse_res.unwrap()
    assert len(parsed.gene_sets) == len(coll.gene_sets)
    assert "T_NK" in parsed.gene_sets


def test_signature_scoring(mock_single_cell_adata: ad.AnnData) -> None:
    coll = get_tme_major_lineage_collection()
    z_res = score_geneset_zscore(mock_single_cell_adata, coll)
    assert isinstance(z_res, Success)
    df_z = z_res.unwrap()
    assert "sample_id" in df_z.columns
    assert "T_NK" in df_z.columns
    assert df_z.height == mock_single_cell_adata.n_obs

    auc_res = score_geneset_auc(mock_single_cell_adata, coll)
    assert isinstance(auc_res, Success)
    df_auc = auc_res.unwrap()
    assert "T_NK" in df_auc.columns
    assert df_auc.height == mock_single_cell_adata.n_obs


def test_compute_geneset_overlap() -> None:
    coll = get_tme_major_lineage_collection()
    overlap_res = compute_geneset_overlap(coll, coll)
    assert isinstance(overlap_res, Success)
    df_over = overlap_res.unwrap()
    assert "jaccard_similarity" in df_over.columns
    # Diagonal self-overlap should be 1.0
    diag = df_over.filter(df_over["signature_a"] == df_over["signature_b"])
    for val in diag["jaccard_similarity"]:
        assert val == 1.0
