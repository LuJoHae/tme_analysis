"""Dataset summary generator and Altair SVG comparison chart creator for ICB single-cell studies."""

from pathlib import Path
import altair as alt
import polars as pl
from returns.result import Result, Success, Failure
from .config import DatasetSpec, DataDirectories, TIER_1_DATASETS


# All 27 datasets evaluated in ICB single-cell & spatial transcriptomics papers.md
ALL_EVALUATED_DATASETS = pl.DataFrame([
    {
        "id": 1, "first_author": "Sade-Feldman", "year": 2018, "journal": "Cell",
        "cancer_type": "Melanoma", "therapy": "Anti-PD-1, anti-CTLA-4", "modality": "scRNA-seq (Smart-seq2)",
        "cell_isolation": "FACS Filtered (CD45+ Immune)",
        "public_status": "Yes", "scrna_seq": "Yes", "objective_response": "Yes (RECIST CR/PR vs SD/PD)",
        "tier": "Tier 1 (Benchmark)", "accession": "GEO GSE120575", "approx_cells": 16291, "approx_patients": 32
    },
    {
        "id": 2, "first_author": "Jerby-Arnon", "year": 2018, "journal": "Cell",
        "cancer_type": "Melanoma", "therapy": "Anti-PD-1, anti-CTLA-4", "modality": "scRNA-seq (Smart-seq2)",
        "cell_isolation": "FACS Filtered (Malignant, CD45+, Stroma)",
        "public_status": "Yes", "scrna_seq": "Yes", "objective_response": "Yes (Prior ICB response)",
        "tier": "Tier 1 (Eligible)", "accession": "GEO GSE115978", "approx_cells": 7186, "approx_patients": 31
    },
    {
        "id": 3, "first_author": "Li", "year": 2019, "journal": "Cell",
        "cancer_type": "Melanoma", "therapy": "ICI", "modality": "scRNA-seq (Smart-seq2)",
        "cell_isolation": "FACS Filtered (CD3+ / CD8+ T cells)",
        "public_status": "Yes", "scrna_seq": "Yes", "objective_response": "Yes (Response documented)",
        "tier": "Tier 1 (Eligible)", "accession": "GEO GSE123139", "approx_cells": 4645, "approx_patients": 25
    },
    {
        "id": 4, "first_author": "de Andrade", "year": 2019, "journal": "JCI Insight",
        "cancer_type": "Melanoma", "therapy": "ICI", "modality": "scRNA-seq",
        "cell_isolation": "FACS Filtered (CD45+ TILs)",
        "public_status": "Yes", "scrna_seq": "Yes", "objective_response": "No (Single Class: PD only)",
        "tier": "Tier 4 (Ineligible)", "accession": "GEO GSE139249", "approx_cells": 2100, "approx_patients": 5
    },
    {
        "id": 5, "first_author": "Pozniak", "year": 2024, "journal": "Cell",
        "cancer_type": "Melanoma", "therapy": "ICB", "modality": "scRNA-seq + Spatial (10X)",
        "cell_isolation": "All Cells (Unsorted Total Tumor Digest)",
        "public_status": "Yes", "scrna_seq": "Yes", "objective_response": "Yes (Responders vs NR)",
        "tier": "Tier 1 (Eligible)", "accession": "KU Leuven RDR", "approx_cells": 45000, "approx_patients": 28
    },
    {
        "id": 6, "first_author": "Alvarez-Breckenridge", "year": 2022, "journal": "Cancer Immunol Res",
        "cancer_type": "Melanoma (Brain Mets)", "therapy": "ICI", "modality": "scRNA-seq + scTCR-seq",
        "cell_isolation": "FACS Filtered (CD45+ & Malignant)",
        "public_status": "Yes", "scrna_seq": "Yes", "objective_response": "Yes (Responders vs NR)",
        "tier": "Tier 1 (Eligible)", "accession": "Broad SCP1493", "approx_cells": 12500, "approx_patients": 24
    },
    {
        "id": 7, "first_author": "Yost", "year": 2019, "journal": "Nature Medicine",
        "cancer_type": "BCC / SCC", "therapy": "Anti-PD-1", "modality": "scRNA-seq + scTCR-seq (10X)",
        "cell_isolation": "FACS Filtered (CD45+ Immune)",
        "public_status": "Yes", "scrna_seq": "Yes", "objective_response": "Yes (Tumor regression / R vs NR)",
        "tier": "Tier 1 (Eligible)", "accession": "GSE123813 / GSE123814", "approx_cells": 79046, "approx_patients": 32
    },
    {
        "id": 8, "first_author": "Bassez", "year": 2021, "journal": "Nature Medicine",
        "cancer_type": "Breast (TNBC/HER2+/ER+)", "therapy": "Anti-PD-1", "modality": "scRNA-seq + CITE-seq (10X)",
        "cell_isolation": "FACS Filtered (Live CD45+ Immune)",
        "public_status": "Partial", "scrna_seq": "Yes", "objective_response": "Surrogate (T-cell expansion)",
        "tier": "Tier 2 (Eligible)", "accession": "Lambrechts / EGA", "approx_cells": 80000, "approx_patients": 40
    },
    {
        "id": 9, "first_author": "Zhang", "year": 2021, "journal": "Cancer Cell",
        "cancer_type": "TNBC", "therapy": "Atezolizumab + chemo", "modality": "scRNA-seq (10X)",
        "cell_isolation": "FACS Filtered (CD45+ & CD45- Fractions)",
        "public_status": "Yes", "scrna_seq": "Yes", "objective_response": "Yes (RECIST CR/PR vs SD/PD)",
        "tier": "Tier 1 (Eligible)", "accession": "GEO GSE169246", "approx_cells": 55000, "approx_patients": 22
    },
    {
        "id": 10, "first_author": "Bi", "year": 2021, "journal": "Cancer Cell",
        "cancer_type": "ccRCC", "therapy": "ICI", "modality": "scRNA-seq (10X)",
        "cell_isolation": "All Cells (Unsorted Total + CD45+ Enriched)",
        "public_status": "Yes", "scrna_seq": "Yes", "objective_response": "Yes (PR vs SD)",
        "tier": "Tier 1 (Eligible)", "accession": "Broad SCP1288", "approx_cells": 34000, "approx_patients": 13
    },
    {
        "id": 11, "first_author": "Ma", "year": 2019, "journal": "Cancer Cell",
        "cancer_type": "HCC / iCCA", "therapy": "ICI / Immunotherapy", "modality": "scRNA-seq (10X)",
        "cell_isolation": "FACS Filtered (CD45+ & CD45- Sorted)",
        "public_status": "Yes", "scrna_seq": "Yes", "objective_response": "Yes (RECIST response)",
        "tier": "Tier 1 (Eligible)", "accession": "GEO GSE125449", "approx_cells": 28000, "approx_patients": 33
    },
    {
        "id": 12, "first_author": "Liu", "year": 2022, "journal": "Nature Cancer",
        "cancer_type": "NSCLC", "therapy": "Anti-PD-1 + chemo", "modality": "scRNA-seq + scTCR-seq (10X)",
        "cell_isolation": "FACS Filtered (CD45+ T cells)",
        "public_status": "Yes", "scrna_seq": "Yes", "objective_response": "Yes (Responsive vs NR)",
        "tier": "Tier 1 (Eligible)", "accession": "GEO GSE179994", "approx_cells": 62000, "approx_patients": 36
    },
    {
        "id": 13, "first_author": "Liu", "year": 2025, "journal": "Cell",
        "cancer_type": "NSCLC", "therapy": "Neoadjuvant Anti-PD-1 + chemo", "modality": "scRNA-seq + scTCR-seq",
        "cell_isolation": "All Cells (Unsorted Total Digest + CD45+ Enriched)",
        "public_status": "Yes", "scrna_seq": "Yes", "objective_response": "Yes (Neoadjuvant MPR/pCR)",
        "tier": "Tier 2 (Eligible)", "accession": "GEO GSE243013", "approx_cells": 120000, "approx_patients": 234
    },
    {
        "id": 14, "first_author": "Luoma", "year": 2022, "journal": "Cell",
        "cancer_type": "HNSCC", "therapy": "Neoadjuvant Anti-PD-1 ± CTLA-4", "modality": "scRNA-seq + scTCR-seq",
        "cell_isolation": "FACS Filtered (CD45+ Immune & Blood T cells)",
        "public_status": "Yes", "scrna_seq": "Yes", "objective_response": "Yes (Neoadjuvant regression)",
        "tier": "Tier 2 (Eligible)", "accession": "GEO GSE200996", "approx_cells": 40000, "approx_patients": 29
    },
    {
        "id": 15, "first_author": "Li", "year": 2023, "journal": "Cancer Cell",
        "cancer_type": "CRC (dMMR/MSI-H)", "therapy": "Neoadjuvant Anti-PD-1", "modality": "scRNA-seq",
        "cell_isolation": "FACS Filtered (CD45+ Immune)",
        "public_status": "Yes", "scrna_seq": "Yes", "objective_response": "Yes (Neoadjuvant pCR)",
        "tier": "Tier 2 (Eligible)", "accession": "GEO GSE205506", "approx_cells": 30000, "approx_patients": 19
    },
    {
        "id": 16, "first_author": "Wu", "year": 2023, "journal": "BMC Medicine",
        "cancer_type": "CRC (MSI-H)", "therapy": "Anti-PD-1", "modality": "scRNA-seq (10X)",
        "cell_isolation": "FACS Filtered (CD45+ Immune)",
        "public_status": "Yes", "scrna_seq": "Yes", "objective_response": "Yes (Resistant vs Sensitive)",
        "tier": "Tier 1 (Eligible)", "accession": "NCBI PRJNA932556", "approx_cells": 38000, "approx_patients": 16
    },
    {
        "id": 17, "first_author": "Hwang", "year": 2022, "journal": "Nature Genetics",
        "cancer_type": "PDAC", "therapy": "Nivolumab (subset N=7)", "modality": "snRNA-seq + Spatial",
        "cell_isolation": "All Nuclei (Unsorted snRNA-seq)",
        "public_status": "Yes", "scrna_seq": "Yes", "objective_response": "Partial",
        "tier": "Tier 2 (Eligible)", "accession": "GEO GSE202051", "approx_cells": 224988, "approx_patients": 7
    },
    {
        "id": 18, "first_author": "Meylan", "year": 2022, "journal": "Immunity",
        "cancer_type": "ccRCC", "therapy": "Nivolumab ± ipilimumab", "modality": "Spatial Visium spots",
        "cell_isolation": "Spatial Tissue Sections (Unsorted Spots)",
        "public_status": "Yes", "scrna_seq": "No (Spatial only)", "objective_response": "Yes (RECIST response)",
        "tier": "Tier 3 (Spatial)", "accession": "GEO GSE175540", "approx_cells": 0, "approx_patients": 24
    },
    {
        "id": 19, "first_author": "Zhang", "year": 2024, "journal": "Nature Comms",
        "cancer_type": "CRC (dMMR/pMMR)", "therapy": "Neoadjuvant Anti-PD-1", "modality": "Stereo-seq + scRNA-seq",
        "cell_isolation": "Spatial Sections & FACS Filtered",
        "public_status": "Restricted", "scrna_seq": "Yes", "objective_response": "Yes (CR/PR vs SD)",
        "tier": "Tier 4 (Restricted)", "accession": "CNCB GSA PRJCA020107", "approx_cells": 25000, "approx_patients": 15
    },
    {
        "id": 20, "first_author": "Mebane", "year": 2025, "journal": "iScience",
        "cancer_type": "TNBC", "therapy": "Pembrolizumab ± SBRT", "modality": "CosMx + scRNA-seq",
        "cell_isolation": "All Cells (CosMx In Situ Spatial + scRNA)",
        "public_status": "Yes", "scrna_seq": "Yes", "objective_response": "Yes (Clearance)",
        "tier": "Tier 1 (Eligible)", "accession": "Zenodo / GSE246613", "approx_cells": 651683, "approx_patients": 4
    },
    {
        "id": 21, "first_author": "Italiano", "year": 2022, "journal": "Nature Medicine",
        "cancer_type": "Sarcoma", "therapy": "Pembrolizumab + cyclo", "modality": "Spatial GeoMx bulk ROI",
        "cell_isolation": "Spatial Bulk ROI (No Single-Cell)",
        "public_status": "No", "scrna_seq": "No (Spatial ROI)", "objective_response": "Yes (Responders vs PD)",
        "tier": "Tier 4 (Ineligible)", "accession": "Upon request", "approx_cells": 0, "approx_patients": 6
    },
    {
        "id": 22, "first_author": "Larroquette", "year": 2022, "journal": "JITC",
        "cancer_type": "NSCLC", "therapy": "ICI", "modality": "Spatial GeoMx bulk ROI",
        "cell_isolation": "Spatial Bulk ROI (No Single-Cell)",
        "public_status": "No", "scrna_seq": "No (Spatial ROI)", "objective_response": "Yes (Durable benefit)",
        "tier": "Tier 4 (Ineligible)", "accession": "Upon request", "approx_cells": 0, "approx_patients": 16
    },
    {
        "id": 23, "first_author": "Park", "year": 2023, "journal": "Cancer Research",
        "cancer_type": "Gastric", "therapy": "ICI", "modality": "Spatial GeoMx bulk ROI",
        "cell_isolation": "Spatial Bulk ROI (No Single-Cell)",
        "public_status": "No", "scrna_seq": "No (Spatial ROI)", "objective_response": "Yes (5 R vs 7 NR)",
        "tier": "Tier 4 (Ineligible)", "accession": "Abstract only", "approx_cells": 0, "approx_patients": 12
    },
    {
        "id": 24, "first_author": "Peyraud", "year": 2025, "journal": "Cell Reports Med",
        "cancer_type": "NSCLC", "therapy": "ICI", "modality": "Spatial GeoMx bulk ROI",
        "cell_isolation": "Spatial Bulk ROI (No Single-Cell)",
        "public_status": "No", "scrna_seq": "No (Spatial ROI)", "objective_response": "Yes (3 R vs 3 PD)",
        "tier": "Tier 4 (Ineligible)", "accession": "Controlled French Ethics", "approx_cells": 0, "approx_patients": 6
    },
    {
        "id": 25, "first_author": "Liu", "year": 2023, "journal": "J Hepatol",
        "cancer_type": "HCC", "therapy": "Anti-PD-1", "modality": "Spatial Visium spots",
        "cell_isolation": "Spatial Spots (Unsorted Visium)",
        "public_status": "Restricted", "scrna_seq": "No (Spatial Visium)", "objective_response": "Yes (3 R vs 5 NR)",
        "tier": "Tier 3 (Spatial)", "accession": "CNCB GSA-Human", "approx_cells": 0, "approx_patients": 8
    },
    {
        "id": 26, "first_author": "Krishna", "year": 2021, "journal": "Cancer Cell",
        "cancer_type": "ccRCC", "therapy": "Nivolumab ± ipilimumab", "modality": "scRNA-seq + Spatial",
        "cell_isolation": "FACS Filtered (CD45+ & CD45- Fractions)",
        "public_status": "Yes", "scrna_seq": "Yes", "objective_response": "Yes (Therapy response)",
        "tier": "Tier 1 (Eligible)", "accession": "PRJNA705464", "approx_cells": 52000, "approx_patients": 21
    },
    {
        "id": 27, "first_author": "Gondal", "year": 2025, "journal": "Scientific Data",
        "cancer_type": "9 Cancer Types", "therapy": "ICB (Anti-PD-1, CTLA-4)", "modality": "Integrated Resource",
        "cell_isolation": "Multi-Cohort (FACS Filtered & Unsorted)",
        "public_status": "Yes", "scrna_seq": "Yes", "objective_response": "Yes (Harmonized R vs NR)",
        "tier": "Tier 1 (Meta-Resource)", "accession": "Zenodo / CELLxGENE", "approx_cells": 355941, "approx_patients": 223
    }
])


def generate_dataset_plots(df: pl.DataFrame, out_dir: Path) -> Result[tuple[Path, Path], str]:
    """Generates Altair charts and exports them as SVG files."""
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        py_df = df.to_pandas()

        tier1_df = py_df[py_df["tier"].str.startswith("Tier 1")]
        
        chart_patients = (
            alt.Chart(tier1_df)
            .mark_bar(color="#2b5c8f")
            .encode(
                x=alt.X("approx_patients:Q", title="Number of Patients"),
                y=alt.Y("first_author:N", sort="-x", title="Dataset / First Author"),
                color=alt.Color("cell_isolation:N", title="Cell Isolation Strategy"),
                tooltip=["first_author", "cancer_type", "cell_isolation", "approx_patients", "approx_cells"]
            )
            .properties(
                title="Tier 1 ICB Single-Cell Datasets: Patient Scale & Cell Isolation",
                width=650,
                height=350
            )
        )

        chart_tiers = (
            alt.Chart(py_df)
            .mark_bar()
            .encode(
                x=alt.X("count():Q", title="Number of Studies"),
                y=alt.Y("tier:N", title="Eligibility Tier"),
                color=alt.Color("tier:N", legend=None),
                tooltip=["tier", "count()"]
            )
            .properties(
                title="Systematic Multi-Criteria Tier Evaluation Summary",
                width=600,
                height=250
            )
        )

        patients_svg = out_dir / "tier1_patient_counts.svg"
        tiers_svg = out_dir / "dataset_tier_summary.svg"

        chart_patients.save(str(patients_svg))
        chart_tiers.save(str(tiers_svg))

        return Success((patients_svg, tiers_svg))
    except Exception as e:
        return Failure(f"Failed to generate Altair SVG plots: {str(e)}")


def build_markdown_tables(df: pl.DataFrame) -> str:
    """Generates markdown comparison tables from Polars dataframe."""
    lines = [
        "## Systematic Multi-Criteria ICB Dataset Evaluation & Cell Isolation Strategy",
        "",
        "| # | First Author | Year | Journal | Cancer Type | Cell Isolation / Pre-filtering | All Cells vs Sorted | Public & Open? | Eligibility Tier | Data Accession |",
        "|:---:|:---|:---:|:---|:---|:---|:---:|:---:|:---|:---|"
    ]
    for row in df.iter_rows(named=True):
        isolation_type = "Unsorted Total" if "All Cells" in row["cell_isolation"] or "Spatial" in row["cell_isolation"] else "FACS Sorted"
        lines.append(
            f"| {row['id']} | {row['first_author']} | {row['year']} | {row['journal']} | "
            f"{row['cancer_type']} | {row['cell_isolation']} | **{isolation_type}** | "
            f"{row['public_status']} | **{row['tier']}** | {row['accession']} |"
        )
    return "\n".join(lines)


def generate_summary_report(dirs: DataDirectories) -> Result[tuple[pl.DataFrame, Path, tuple[Path, Path]], str]:
    """Generates comprehensive Polars summary dataframe, markdown tables, and SVG plots."""
    dirs.reports_dir.mkdir(parents=True, exist_ok=True)
    report_md_path = dirs.reports_dir / "icb_single_cell_datasets_summary.md"

    match generate_dataset_plots(ALL_EVALUATED_DATASETS, dirs.reports_dir):
        case Success(plot_paths):
            try:
                md_tables = build_markdown_tables(ALL_EVALUATED_DATASETS)
                report_content = f"# ICB Single-Cell & Spatial Transcriptomics Dataset Evaluation Report\n\n{md_tables}\n"
                report_md_path.write_text(report_content, encoding="utf-8")
                return Success((ALL_EVALUATED_DATASETS, report_md_path, plot_paths))
            except Exception as e:
                return Failure(f"Failed to write summary report markdown: {str(e)}")
        case Failure(err):
            return Failure(err)
