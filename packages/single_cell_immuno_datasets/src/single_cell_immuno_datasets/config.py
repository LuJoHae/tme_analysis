"""Configuration specifications for ICB single-cell and spatial transcriptomics datasets."""

from pathlib import Path
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Some, Nothing


class DatasetSpec(BaseModel):
    """Immutable metadata specification for an ICB single-cell dataset."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    accession: str
    paper_title: str
    first_author: str
    year: int
    journal: str
    doi: str
    cancer_type: str
    therapy: str
    modality: str
    sequencing_tech: str
    tier: str
    public_status: str
    objective_response: str
    download_urls: dict[str, str]
    obs_mapping: dict[str, str]
    # Hardcoded patient_id → response label (used when response metadata is not downloadable)
    response_map: dict[str, str] = {}
    cell_count_approx: Maybe[int] = Nothing
    patient_count_approx: Maybe[int] = Nothing


class QualityControlSpec(BaseModel):
    """Quality control threshold specification for single-cell preprocessing."""
    model_config = ConfigDict(frozen=True)
    
    min_genes: int = 300
    max_genes: int = 6000
    min_counts: int = 500
    max_counts: int = 30000
    max_mt_content: float = 20.0


class DataDirectories(BaseModel):
    """Immutable target directory configuration defaulting to remote server storage."""
    model_config = ConfigDict(frozen=True)
    
    base_dir: Path = Path("/storage/halu/data-test")
    raw_dir: Path = Path("/storage/halu/data-test/raw")
    preprocessed_dir: Path = Path("/storage/halu/data-test/preprocessed")
    reports_dir: Path = Path("/storage/halu/data-test/reports")


    @classmethod
    def with_base(cls, custom_base: Path) -> "DataDirectories":
        return cls(
            base_dir=custom_base,
            raw_dir=custom_base / "raw",
            preprocessed_dir=custom_base / "preprocessed",
            reports_dir=custom_base / "reports",
        )


# Registry of Tier 1 (Benchmark & Primary Eligible ICB scRNA-seq Datasets)
TIER_1_DATASETS: tuple[DatasetSpec, ...] = (
    DatasetSpec(
        accession="GSE120575",
        paper_title="Defining T cell states associated with response to checkpoint immunotherapy in melanoma",
        first_author="Sade-Feldman",
        year=2018,
        journal="Cell",
        doi="10.1016/j.cell.2018.10.038",
        cancer_type="Melanoma",
        therapy="Anti-PD-1, anti-CTLA-4",
        modality="scRNA-seq",
        sequencing_tech="Smart-seq2",
        tier="Tier 1 (Benchmark)",
        public_status="Yes",
        objective_response="Yes (RECIST CR/PR vs SD/PD)",
        download_urls={
            "tpm": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz",
            "patient_meta": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_patient_ID_single_cells.txt.gz",
        },
        obs_mapping={"Cell_ID": "original.barcode", "Patient_ID": "patient"},
        cell_count_approx=Some(16291),
        patient_count_approx=Some(32),
    ),
    DatasetSpec(
        accession="GSE115978",
        paper_title="A melanoma cell state associated with immune exclusion and resistance to immunotherapy",
        first_author="Jerby-Arnon",
        year=2018,
        journal="Cell",
        doi="10.1016/j.cell.2018.09.006",
        cancer_type="Melanoma",
        therapy="Anti-PD-1, anti-CTLA-4",
        modality="scRNA-seq",
        sequencing_tech="Smart-seq2",
        tier="Tier 1 (Eligible)",
        public_status="Yes",
        objective_response="Yes (Prior ICB response)",
        download_urls={
            "tpm": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115978/suppl/GSE115978_tpm.csv.gz",
            "meta": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115978/suppl/GSE115978_cell.annotations.csv.gz",
        },
        obs_mapping={"samples": "sample", "Cohort": "cohort", "cell.types": "cell_type"},
        cell_count_approx=Some(7186),
        patient_count_approx=Some(31),
    ),
    DatasetSpec(
        accession="GSE123139",
        paper_title="Dysfunctional CD8 T cells form a proliferative pool supporting response to immune checkpoint blockade",
        first_author="Li",
        year=2019,
        journal="Cell",
        doi="10.1016/j.cell.2019.08.004",
        cancer_type="Melanoma",
        therapy="ICI",
        modality="scRNA-seq",
        sequencing_tech="Smart-seq2",
        tier="Tier 1 (Eligible)",
        public_status="Yes",
        objective_response="Yes (Response documented)",
        download_urls={
            "raw_tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123139/suppl/GSE123139_RAW.tar",
        },
        obs_mapping={},
        cell_count_approx=Some(4645),
        patient_count_approx=Some(25),
    ),
    DatasetSpec(
        accession="GSE123813",
        paper_title="Clonal replacement of tumor-specific T cells following PD-1 blockade",
        first_author="Yost",
        year=2019,
        journal="Nature Medicine",
        doi="10.1038/s41591-019-0522-3",
        cancer_type="BCC / SCC",
        therapy="Anti-PD-1",
        modality="scRNA-seq + scTCR-seq",
        sequencing_tech="10x_UMI",
        tier="Tier 1 (Eligible)",
        public_status="Yes",
        objective_response="Yes (Tumor regression / R vs NR)",
        download_urls={
            "bcc_counts": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123813/suppl/GSE123813_bcc_scRNA_counts.txt.gz",
            "scc_counts": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123813/suppl/GSE123813_scc_scRNA_counts.txt.gz",
        },
        obs_mapping={},
        response_map={
            "su001": "Responder",
            "su002": "Responder",
            "su003": "Responder",
            "su004": "Responder",
            "su005": "Non-responder",
            "su006": "Non-responder",
            "su007": "Non-responder",
            "su008": "Non-responder",
            "su009": "Responder",
            "su010": "Non-responder",
            "su011": "Responder",
            "su012": "Responder",
            "su013": "Non-responder",
            "su014": "Non-responder",
        },
        cell_count_approx=Some(79046),
        patient_count_approx=Some(32),
    ),
    DatasetSpec(
        accession="GSE169246",
        paper_title="Single-cell analyses reveal key immune cell subsets associated with response to PD-L1 blockade in TNBC",
        first_author="Zhang",
        year=2021,
        journal="Cancer Cell",
        doi="10.1016/j.ccell.2021.09.010",
        cancer_type="TNBC",
        therapy="Atezolizumab + paclitaxel",
        modality="scRNA-seq",
        sequencing_tech="10x_UMI",
        tier="Tier 1 (Eligible)",
        public_status="Yes",
        objective_response="Yes (RECIST CR/PR vs SD/PD)",
        download_urls={
            "counts": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169246/suppl/GSE169246_TNBC_RNA.counts.mtx.gz",
            "barcodes": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169246/suppl/GSE169246_TNBC_RNA.barcode.tsv.gz",
            "features": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169246/suppl/GSE169246_TNBC_RNA.feature.tsv.gz",
        },
        obs_mapping={},
        cell_count_approx=Some(55000),
        patient_count_approx=Some(22),
    ),
    DatasetSpec(
        accession="GSE159115",
        paper_title="Tumor and immune cell dynamics in clear cell renal cell carcinoma during immune checkpoint blockade",
        first_author="Bi",
        year=2021,
        journal="Cancer Cell",
        doi="10.1016/j.ccell.2021.02.015",
        cancer_type="ccRCC",
        therapy="ICI",
        modality="scRNA-seq",
        sequencing_tech="10x_UMI",
        tier="Tier 1 (Eligible)",
        public_status="Yes",
        objective_response="Yes (PR vs SD)",
        download_urls={
            "raw_tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159115/suppl/GSE159115_RAW.tar",
            "anno": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159115/suppl/GSE159115_ccRCC_anno.csv.gz",
        },
        obs_mapping={},
        # Bi 2021 Cancer Cell: patient column in h5ad = GSM_SI format from H5 filenames.
        # Anno.csv provides: SI_18854=SS_2005(PR), SI_18856=SS_2005-post(PR),
        # SI_22604=SS_2022(PR), SI_22605=SS_2022-post(PR),
        # SI_23459=SS_2023(SD), SI_23843=SS_2026(SD).
        # Post-treatment samples inherit the same response as the pre-treatment from same patient.
        response_map={
            "GSM4819725_SI_18854": "PR",   # SS_2005 pre
            "GSM4819726_SI_18856": "PR",   # SS_2005 post
            "GSM4819736_SI_22604": "PR",   # SS_2022 pre
            "GSM4819735_SI_22605": "PR",   # SS_2022 post
            "GSM4819737_SI_23459": "SD",   # SS_2023 pre
            "GSM4819738_SI_23843": "SD",   # SS_2026 pre (no matching post in anno)
            "GSM4819734_SI_22368": "SD",   # SS_2017 pre (anno confirms mapping)
            "GSM4819733_SI_22369": "SD",   # SS_2017 post
        },
        cell_count_approx=Some(34000),
        patient_count_approx=Some(13),
    ),

    DatasetSpec(
        accession="GSE125449",
        paper_title="Single-cell transcriptomic analysis defines heterogeneity and immune cell states in liver cancers",
        first_author="Ma",
        year=2019,
        journal="Cancer Cell",
        doi="10.1016/j.ccell.2019.08.007",
        cancer_type="HCC / iCCA",
        therapy="ICI / Immunotherapy",
        modality="scRNA-seq",
        sequencing_tech="10x_UMI",
        tier="Tier 1 (Eligible)",
        public_status="Yes",
        objective_response="Yes (RECIST response)",
        download_urls={
            "set1_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125449/suppl/GSE125449_Set1_matrix.mtx.gz",
            "set1_barcodes": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125449/suppl/GSE125449_Set1_barcodes.tsv.gz",
            "set1_genes": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125449/suppl/GSE125449_Set1_genes.tsv.gz",
        },
        obs_mapping={},
        cell_count_approx=Some(28000),
        patient_count_approx=Some(33),
    ),
    DatasetSpec(
        accession="GSE179994",
        paper_title="Precursor exhausted T cells expand upon anti-PD-1 blockade in non-small cell lung cancer",
        first_author="Liu",
        year=2022,
        journal="Nature Cancer",
        doi="10.1038/s43018-021-00292-8",
        cancer_type="NSCLC",
        therapy="Anti-PD-1 + chemo",
        modality="scRNA-seq + scTCR-seq",
        sequencing_tech="10x_UMI",
        tier="Tier 1 (Eligible)",
        public_status="Yes",
        objective_response="Yes (Responsive vs NR)",
        download_urls={
            "counts": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179994/suppl/GSE179994_all.Tcell.rawCounts.rds.gz",
            "meta": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179994/suppl/GSE179994_Tcell.metadata.tsv.gz",
        },
        obs_mapping={},
        cell_count_approx=Some(62000),
        patient_count_approx=Some(36),
    ),
    DatasetSpec(
        accession="PRJNA932556",
        paper_title="IL-1β-driven MDSC infiltration suppresses CD8 T cells and enhances resistance to anti-PD-1 in mCRC",
        first_author="Wu",
        year=2023,
        journal="BMC Medicine",
        doi="10.1186/s12916-023-02866-y",
        cancer_type="CRC (MSI-H)",
        therapy="Anti-PD-1",
        modality="scRNA-seq",
        sequencing_tech="10x_UMI",
        tier="Tier 1 (Eligible)",
        public_status="Yes (NCBI SRA)",
        objective_response="Yes (Resistant vs Sensitive)",
        download_urls={
            "sra": "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA932556",
        },
        obs_mapping={},
        cell_count_approx=Some(38000),
        patient_count_approx=Some(16),
    ),
    DatasetSpec(
        accession="GSE171306",
        paper_title="Single-cell transcriptomics reveals multi-regional immune landscapes in ccRCC",
        first_author="Krishna",
        year=2021,
        journal="Cancer Cell",
        doi="10.1016/j.ccell.2021.03.007",
        cancer_type="ccRCC",
        therapy="Nivolumab ± ipilimumab",
        modality="scRNA-seq + Spatial",
        sequencing_tech="10x_UMI",
        tier="Tier 1 (Eligible)",
        public_status="Yes",
        objective_response="Yes (Therapy efficacy)",
        download_urls={
            "raw_tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE171nnn/GSE171306/suppl/GSE171306_RAW.tar",
        },
        obs_mapping={},
        # Krishna 2021 Cancer Cell: GEO supplement contains ccRCC1 (responder) and ccRCC2 (non-responder)
        response_map={
            "GSM5222644_ccRCC1": "Responder",
            "GSM5222645_ccRCC2": "Non-responder",
        },
        cell_count_approx=Some(52000),
        patient_count_approx=Some(21),
    ),

    DatasetSpec(
        accession="Gondal2025",
        paper_title="An integrated single-cell RNA-seq resource of immune checkpoint blockade-treated cancer patients",
        first_author="Gondal",
        year=2025,
        journal="Scientific Data",
        doi="10.1038/s41597-025-04381-6",
        cancer_type="9 Cancer Types (Pan-Cancer)",
        therapy="ICB (Anti-PD-1, CTLA-4)",
        modality="scRNA-seq Integrated Resource",
        sequencing_tech="Multi-platform",
        tier="Tier 1 (Meta-Resource)",
        public_status="Yes (Zenodo / CELLxGENE)",
        objective_response="Yes (Harmonized R vs NR)",
        download_urls={
            "zenodo": "https://zenodo.org/api/records/10407126/files/Rshinydata_singlecell-20231219T155916Z-001.zip/content",
        },
        obs_mapping={},
        cell_count_approx=Some(355941),
        patient_count_approx=Some(223),
    ),
)
