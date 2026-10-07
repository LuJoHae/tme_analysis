"""Unit tests for the 22 Tier 0 Premier Benchmark Core scRNA-seq cohort loaders and preprocessing."""

from pathlib import Path
import tempfile
import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from returns.maybe import Some
from returns.result import Failure, Success

from tme_datasets.models import QualityControlSpec
from tme_datasets.preprocessing.normalization import standardize_dual_layers, compute_log1p_norm_matrix
from tme_datasets.preprocessing.qc import apply_quality_control
from tme_datasets.providers.tier0_single_cell import (
    load_gse246613_breast,
    load_gse300475_breast,
    load_gse212707_breast,
    load_gse236581_crc,
    load_gse299651_crc,
    load_cellxgene_829a3cd1_crc,
    load_gse270680_gastric,
    load_gse313642_hcc,
    load_gse245906_hcc,
    load_gse301741_hnscc,
    load_gse287301_hnscc,
    load_gse200996_hnscc,
    load_cellxgene_7b20c613_melanoma,
    load_cellxgene_05a8c945_crc,
    load_cellxgene_6f9de485_breast,
    load_gse218429_melanoma,
    load_gse344166_melanoma,
    load_gse207422_nsclc,
    load_gse317309_nsclc,
    load_gse243013_nsclc,
    load_gse233203_nsclc,
    load_gse311789_pdac,
    load_gse316195_pdac,
    load_gse210038_ccrcc,
    load_gse314072_ccrcc,
    _harmonize_single_response,
)
from tme_datasets.registry import (
    DATASET_ALIASES,
    DATASET_REGISTRY,
    get_dataset_spec,
    resolve_dataset_id,
)

TIER0_DATASET_IDS = [
    "GSE246613",
    "GSE300475",
    "GSE212707",
    "GSE236581",
    "GSE299651",
    "CELLxGENE_829a3cd1",
    "GSE270680",
    "GSE313642",
    "GSE245906",
    "GSE301741",
    "GSE287301",
    "GSE200996",
    "CELLxGENE_7b20c613",
    "CELLxGENE_05a8c945",
    "CELLxGENE_6f9de485",
    "GSE218429",
    "GSE344166",
    "GSE207422",
    "GSE317309",
    "GSE243013",
    "GSE233203",
    "GSE311789",
    "GSE316195",
    "GSE210038",
    "GSE314072",
]


def test_tier0_registry_coverage() -> None:
    """Verify all 25 Tier 0 benchmark datasets are properly registered with specs."""
    for dataset_id in TIER0_DATASET_IDS:
        assert dataset_id in DATASET_REGISTRY, f"{dataset_id} missing from DATASET_REGISTRY"
        spec_maybe = get_dataset_spec(dataset_id)
        assert isinstance(spec_maybe, Some)
        spec = spec_maybe.unwrap()
        assert spec.id == dataset_id
        assert spec.cancer_type != ""
        assert spec.has_response_labels is True
        assert isinstance(spec.qc_spec, Some), f"{dataset_id} missing QualityControlSpec"

    # snRNA-seq GSE316195 has strict 10% MT cutoff
    pdac_sn_spec = get_dataset_spec("GSE316195").unwrap()
    qc = pdac_sn_spec.qc_spec.unwrap()
    assert qc.max_pct_mitochondrial == 10.0


def test_tier0_alias_resolution() -> None:
    """Verify aliases resolve accurately to canonical Tier 0 identifiers."""
    assert resolve_dataset_id("GSE246613-Breast") == "GSE246613"
    assert resolve_dataset_id("GSE236581-CRC") == "GSE236581"
    assert resolve_dataset_id("GSE200996-HNSCC") == "GSE200996"
    assert resolve_dataset_id("GSE218429-Melanoma") == "GSE218429"
    assert resolve_dataset_id("GSE316195-PDAC") == "GSE316195"
    assert resolve_dataset_id("CELLxGENE-829a3cd1") == "CELLxGENE_829a3cd1"
    assert resolve_dataset_id("CELLxGENE-7b20c613") == "CELLxGENE_7b20c613"
    assert resolve_dataset_id("05a8c945") == "CELLxGENE_05a8c945"
    assert resolve_dataset_id("6f9de485") == "CELLxGENE_6f9de485"
    assert resolve_dataset_id("GSE233203-NSCLC") == "GSE233203"


def test_response_harmonization() -> None:
    """Verify _harmonize_single_response cleanly categorizes diverse clinical response tokens."""
    # Responders
    assert _harmonize_single_response("Favourable")[0] == "responder"
    assert _harmonize_single_response("Complete Response")[0] == "responder"
    assert _harmonize_single_response("CR")[0] == "responder"
    assert _harmonize_single_response("PR")[0] == "responder"
    assert _harmonize_single_response("pCR")[0] == "responder"
    assert _harmonize_single_response("MPR")[0] == "responder"
    assert _harmonize_single_response("Response")[0] == "responder"
    assert _harmonize_single_response("High")[0] == "responder"

    # Non-responders
    assert _harmonize_single_response("Unfavourable")[0] == "non-responder"
    assert _harmonize_single_response("Progressive Disease")[0] == "non-responder"
    assert _harmonize_single_response("PD")[0] == "non-responder"
    assert _harmonize_single_response("RD")[0] == "non-responder"
    assert _harmonize_single_response("non-MPR")[0] == "non-responder"
    assert _harmonize_single_response("NMPR")[0] == "non-responder"
    assert _harmonize_single_response("Non-response")[0] == "non-responder"
    assert _harmonize_single_response("Low")[0] == "non-responder"

    # Stable
    assert _harmonize_single_response("Stable Disease")[0] == "stable"
    assert _harmonize_single_response("SD")[0] == "stable"
    assert _harmonize_single_response("Medium")[0] == "stable"

    # Not evaluable / Unknown
    assert _harmonize_single_response("NE")[0] == "not-evaluable"
    assert _harmonize_single_response("NA")[0] == "not-evaluable"
    assert _harmonize_single_response(None)[0] == "not-evaluable"


def test_tier0_loaders_missing_raw_files() -> None:
    """Verify loaders fail gracefully with informative error messages when raw files are absent."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        fake_path = Path(tmp_dir)

        loaders = [
            load_gse246613_breast,
            load_gse300475_breast,
            load_gse212707_breast,
            load_gse236581_crc,
            load_gse299651_crc,
            load_cellxgene_829a3cd1_crc,
            load_gse270680_gastric,
            load_gse313642_hcc,
            load_gse245906_hcc,
            load_gse301741_hnscc,
            load_gse287301_hnscc,
            load_gse200996_hnscc,
            load_cellxgene_7b20c613_melanoma,
            load_cellxgene_05a8c945_crc,
            load_cellxgene_6f9de485_breast,
            load_gse218429_melanoma,
            load_gse344166_melanoma,
            load_gse207422_nsclc,
            load_gse317309_nsclc,
            load_gse243013_nsclc,
            load_gse233203_nsclc,
            load_gse311789_pdac,
            load_gse316195_pdac,
            load_gse210038_ccrcc,
            load_gse314072_ccrcc,
        ]

        for loader in loaders:
            res = loader(fake_path)
            assert isinstance(res, Failure), f"{loader.__name__} did not fail on missing data"
            assert "not found" in res.failure().lower() or "missing" in res.failure().lower()


def test_quality_control_filtering() -> None:
    """Test QualityControl filtering on low-quality cells, dead droplets, and high-mitochondrial noise."""
    # Create synthetic AnnData with 5 cells:
    # Cell 0: Good cell (counts=1000, 300 genes, MT=5%)
    # Cell 1: Low count cell (counts=300 < 500) -> should be filtered
    # Cell 2: Low gene count (counts=1000, 100 genes < 200) -> should be filtered
    # Cell 3: High MT cell (counts=1000, 300 genes, MT=40% > 20%) -> should be filtered
    # Cell 4: Good cell (counts=1200, 250 genes, MT=8%)
    n_genes = 400
    X = sp.lil_matrix((5, n_genes), dtype=np.float32)

    # Cell 0
    X[0, :300] = 3
    X[0, 390:400] = 10  # 100 MT counts out of 1000 = 10%

    # Cell 1: total counts 100
    X[1, :50] = 2

    # Cell 2: 100 genes with large values (counts=1000, but only 100 genes)
    X[2, :100] = 10

    # Cell 3: 40% MT counts
    X[3, :200] = 3  # 600 counts
    X[3, 380:400] = 20  # 400 counts MT -> 400/1000 = 40%

    # Cell 4
    X[4, :250] = 4
    X[4, 395:400] = 10  # 50 counts MT out of 1050 ~ 4.7%

    var_names = [f"GENE_{i}" for i in range(380)] + [f"MT-ND{i}" for i in range(20)]
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(5)])
    var = pd.DataFrame(index=var_names)

    adata = ad.AnnData(X=X.tocsr(), obs=obs, var=var)

    qc_spec = QualityControlSpec(
        min_genes_per_cell=200,
        max_genes_per_cell=8000,
        min_counts_per_cell=500,
        max_pct_mitochondrial=20.0,
    )

    qc_res = apply_quality_control(adata, qc_spec=qc_spec)
    assert isinstance(qc_res, Success)
    filtered = qc_res.unwrap()

    # Only Cell 0 and Cell 4 should pass
    assert filtered.n_obs == 2
    assert list(filtered.obs_names) == ["cell_0", "cell_4"]
    assert "total_counts" in filtered.obs.columns
    assert "n_genes_by_counts" in filtered.obs.columns
    assert "pct_counts_mt" in filtered.obs.columns
    assert filtered.uns["qc_spec"]["n_cells_pre_qc"] == 5
    assert filtered.uns["qc_spec"]["n_cells_post_qc"] == 2


def test_standardize_dual_layers() -> None:
    """Verify standardize_dual_layers creates accurate float32 sparse CSR layers."""
    # Synthetic AnnData with 3 cells x 4 genes
    raw_X = np.array([
        [10.0, 0.0, 5.0, 0.0],
        [0.0, 20.0, 0.0, 10.0],
        [5.0, 5.0, 5.0, 5.0],
    ], dtype=np.float32)

    adata = ad.AnnData(
        X=raw_X,
        obs=pd.DataFrame(index=["cell_1", "cell_2", "cell_3"]),
        var=pd.DataFrame(index=["GENE1", "GENE2", "GENE3", "GENE4"]),
    )

    res = standardize_dual_layers(adata, target_sum=1e4)
    assert isinstance(res, Success)
    std_adata = res.unwrap()

    # Verify X is sparse CSR float32
    assert sp.isspmatrix_csr(std_adata.X)
    assert std_adata.X.dtype == np.float32
    assert np.allclose(std_adata.X.toarray(), raw_X)

    # Verify layers["counts"] is identical to X
    assert "counts" in std_adata.layers
    assert sp.isspmatrix_csr(std_adata.layers["counts"])
    assert np.allclose(std_adata.layers["counts"].toarray(), raw_X)

    # Verify layers["log1p_norm"]
    assert "log1p_norm" in std_adata.layers
    assert sp.isspmatrix_csr(std_adata.layers["log1p_norm"])
    assert std_adata.layers["log1p_norm"].dtype == np.float32

    # Check cell 0 log1p_norm: total counts = 15, scale = 10000 / 15
    expected_c0_g0 = np.log1p(10.0 * (10000.0 / 15.0))
    actual_c0_g0 = std_adata.layers["log1p_norm"][0, 0]
    assert np.isclose(actual_c0_g0, expected_c0_g0, rtol=1e-5)


def test_mock_loaders_for_8_verified_cohorts() -> None:
    """Verify all 8 verified Tier 0 loaders construct AnnData adhering to the standardized schema."""
    import gzip
    import tarfile
    import scipy.io as sio

    with tempfile.TemporaryDirectory() as tmp_str:
        tmp_dir = Path(tmp_str)

        # 1. Mock CELLxGENE_7b20c613 (Melanoma)
        dir_7b = tmp_dir / "7b"
        dir_7b.mkdir()
        obs_7b = pd.DataFrame({
            "PMID_donor_id": ["P1_d1", "P1_d2"],
            "Combined_outcome": ["Favourable", "Unfavourable"],
            "outcome": ["CR", "PD"],
            "cell_type": ["T cell", "Melanoma"],
        }, index=["c1", "c2"])
        a_7b = ad.AnnData(X=sp.csr_matrix(np.array([[5, 10], [20, 3]], dtype=np.float32)), obs=obs_7b, var=pd.DataFrame(index=["CD8A", "MLANA"]))
        a_7b.write_h5ad(dir_7b / "7b20c613.h5ad")
        res_7b = load_cellxgene_7b20c613_melanoma(dir_7b)
        assert isinstance(res_7b, Success)
        loaded_7b = res_7b.unwrap()
        assert list(loaded_7b.obs["clinical_response"]) == ["responder", "non-responder"]
        assert list(loaded_7b.obs["patient_id"]) == ["P1_d1", "P1_d2"]

        # 2. Mock CELLxGENE_05a8c945 (CRC)
        dir_05 = tmp_dir / "05"
        dir_05.mkdir()
        obs_05 = pd.DataFrame({
            "donor_id": ["D1", "D2"],
            "treatment_response": ["responder", "non-responder"],
            "RECIST": ["CR: complete response", "PD: progressive disease"],
            "treatment_status_before_resection": ["Pre-treatment", "Pre-treatment"],
        }, index=["c1", "c2"])
        a_05 = ad.AnnData(X=sp.csr_matrix(np.array([[12, 1], [3, 14]], dtype=np.float32)), obs=obs_05, var=pd.DataFrame(index=["EPCAM", "CD3D"]))
        a_05.write_h5ad(dir_05 / "05a8c945.h5ad")
        res_05 = load_cellxgene_05a8c945_crc(dir_05)
        assert isinstance(res_05, Success)
        loaded_05 = res_05.unwrap()
        assert list(loaded_05.obs["clinical_response"]) == ["responder", "non-responder"]
        assert list(loaded_05.obs["patient_id"]) == ["D1", "D2"]

        # 3. Mock CELLxGENE_6f9de485 (TNBC)
        dir_6f = tmp_dir / "6f"
        dir_6f.mkdir()
        obs_6f = pd.DataFrame({
            "donor_id": ["TNBC_01", "TNBC_02"],
            "pCR_status": ["pCR", "RD"],
        }, index=["c1", "c2"])
        a_6f = ad.AnnData(X=sp.csr_matrix(np.array([[8, 4], [2, 16]], dtype=np.float32)), obs=obs_6f, var=pd.DataFrame(index=["KRT5", "CD8B"]))
        a_6f.write_h5ad(dir_6f / "6f9de485.h5ad")
        res_6f = load_cellxgene_6f9de485_breast(dir_6f)
        assert isinstance(res_6f, Success)
        loaded_6f = res_6f.unwrap()
        assert list(loaded_6f.obs["clinical_response"]) == ["responder", "non-responder"]
        assert list(loaded_6f.obs["patient_id"]) == ["TNBC_01", "TNBC_02"]

        # 4. Mock GSE207422 (NSCLC UMI + Excel)
        dir_207 = tmp_dir / "207"
        dir_207.mkdir()
        umi_df = pd.DataFrame({
            "Gene": ["CD8A", "TP53"],
            "BD_immune01_c1": [10, 5],
            "BD_immune02_c2": [2, 8],
        })
        umi_path = dir_207 / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
        with gzip.open(umi_path, "wt") as f:
            umi_df.to_csv(f, sep="\t", index=False)
        meta_207 = pd.DataFrame({
            "Sample": ["BD_immune01", "BD_immune02"],
            "Patient": ["P01", "P02"],
            "Pathologic Response": ["MPR", "NMPR"],
            "RECIST": ["PR", "PD"],
            "Chemotherapy": ["Carboplatin", "Carboplatin"],
        })
        meta_207.to_excel(dir_207 / "GSE207422_NSCLC_scRNAseq_metadata.xlsx", index=False)
        res_207 = load_gse207422_nsclc(dir_207)
        assert isinstance(res_207, Success)
        loaded_207 = res_207.unwrap()
        assert list(loaded_207.obs["clinical_response"]) == ["responder", "non-responder"]
        assert list(loaded_207.obs["patient_id"]) == ["P01", "P02"]

        # 5. Mock GSE243013 (NSCLC 10x counts)
        dir_243 = tmp_dir / "243"
        dir_243.mkdir()
        mat_243 = sp.csr_matrix(np.array([[20, 5], [10, 15]], dtype=np.float32))
        sio.mmwrite(str(dir_243 / "GSE243013_NSCLC_immune_scRNA_counts.mtx"), mat_243)
        with open(dir_243 / "GSE243013_NSCLC_immune_scRNA_counts.mtx", "rb") as f_in, gzip.open(dir_243 / "GSE243013_NSCLC_immune_scRNA_counts.mtx.gz", "wb") as f_out:
            f_out.write(f_in.read())
        with gzip.open(dir_243 / "GSE243013_barcodes.csv.gz", "wt") as f:
            f.write("barcode\ncell_1\ncell_2\n")
        with gzip.open(dir_243 / "GSE243013_genes.csv.gz", "wt") as f:
            f.write("geneSymbol\nCD8A\nCD4\n")
        meta_243 = pd.DataFrame({
            "cellID": ["cell_1", "cell_2"],
            "sampleID": ["S1", "S2"],
            "pathological_response": ["MPR", "non-MPR"],
            "major_cell_type": ["T cell", "Myeloid"],
        })
        with gzip.open(dir_243 / "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz", "wt") as f:
            meta_243.to_csv(f, index=False)
        res_243 = load_gse243013_nsclc(dir_243)
        assert isinstance(res_243, Success)
        loaded_243 = res_243.unwrap()
        assert list(loaded_243.obs["clinical_response"]) == ["responder", "non-responder"]
        assert list(loaded_243.obs["sample_id"]) == ["S1", "S2"]

        # Helper to create mock tar with MTX triplet
        def _make_tar_with_mtx(tar_file: Path, prefix: str) -> None:
            sub = tmp_dir / f"sub_{prefix}"
            sub.mkdir(parents=True, exist_ok=True)
            # Matrix: cells x genes
            mat = sp.csr_matrix(np.ones((2, 2), dtype=np.float32) * 50)
            m_path = sub / f"{prefix}_matrix.mtx"
            sio.mmwrite(str(m_path), mat)
            with open(m_path, "rb") as fi, gzip.open(sub / f"{prefix}_matrix.mtx.gz", "wb") as fo:
                fo.write(fi.read())
            with gzip.open(sub / f"{prefix}_barcodes.tsv.gz", "wt") as fo:
                fo.write("AAAC-1\nAAAG-1\n")
            with gzip.open(sub / f"{prefix}_features.tsv.gz", "wt") as fo:
                fo.write("ENSG01\tGENE1\tGene Expression\nENSG02\tGENE2\tGene Expression\n")
            with tarfile.open(tar_file, "w") as tar:
                tar.add(sub / f"{prefix}_matrix.mtx.gz", arcname=f"{prefix}_matrix.mtx.gz")
                tar.add(sub / f"{prefix}_barcodes.tsv.gz", arcname=f"{prefix}_barcodes.tsv.gz")
                tar.add(sub / f"{prefix}_features.tsv.gz", arcname=f"{prefix}_features.tsv.gz")

        # 6. Mock GSE233203 (NSCLC RAW.tar + series matrix)
        dir_233 = tmp_dir / "233"
        dir_233.mkdir()
        _make_tar_with_mtx(dir_233 / "GSE233203_RAW.tar", "GSM7412612_NCCLu_162")
        with gzip.open(dir_233 / "GSE233203_series_matrix.txt.gz", "wt") as f:
            f.write("!Sample_title\t\"NCCLu_162\"\n")
            f.write("!Sample_geo_accession\t\"GSM7412612\"\n")
            f.write("!Sample_characteristics_ch1\t\"therapeutic response: Response\"\n")
        res_233 = load_gse233203_nsclc(dir_233)
        assert isinstance(res_233, Success)
        loaded_233 = res_233.unwrap()
        assert list(loaded_233.obs["clinical_response"]) == ["responder", "responder"]

        # 7. Mock GSE200996 (HNSCC RAW.tar + metadata TSV)
        dir_200 = tmp_dir / "200"
        dir_200.mkdir()
        _make_tar_with_mtx(dir_200 / "GSE200996_RAW.tar", "GSM6048142_P13_post")
        meta_200 = pd.DataFrame({
            "Unnamed: 0": ["AAAC-1", "AAAG-1"],
            "Patient_ID": ["P13", "P13"],
            "Path_response": ["High", "High"],
        })
        with gzip.open(dir_200 / "GSE200996_CD4.tumor.single.cell.meta.data.txt.gz", "wt") as f:
            meta_200.to_csv(f, sep="\t", index=False)
        res_200 = load_gse200996_hnscc(dir_200)
        assert isinstance(res_200, Success)
        loaded_200 = res_200.unwrap()
        assert list(loaded_200.obs["clinical_response"]) == ["responder", "responder"]
        assert list(loaded_200.obs["patient_id"]) == ["P13", "P13"]

        # 8. Mock GSE316195 (PDAC RAW.tar + series matrix)
        dir_316 = tmp_dir / "316"
        dir_316.mkdir()
        _make_tar_with_mtx(dir_316 / "GSE316195_RAW.tar", "GSM9447228_GM18")
        with gzip.open(dir_316 / "GSE316195_series_matrix.txt.gz", "wt") as f:
            f.write("!Sample_title\t\"GM18\"\n")
            f.write("!Sample_geo_accession\t\"GSM9447228\"\n")
            f.write("!Sample_characteristics_ch1\t\"response: PR\"\n")
        res_316 = load_gse316195_pdac(dir_316)
        assert isinstance(res_316, Success)
        loaded_316 = res_316.unwrap()
        assert list(loaded_316.obs["clinical_response"]) == ["responder", "responder"]
        assert loaded_316.uns.get("is_single_nucleus") is True
