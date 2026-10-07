# Comprehensive Catalog: Public Unrestricted Human Solid Tumor scRNA-seq Datasets

**Generated**: 2026-09-23 13:55:37 UTC | **Source**: NCBI GEO & CZ CELLxGENE | **Access**: Public & Unrestricted

> [!NOTE]
> This catalog provides full experimental methodology, sequencing technology/chemistry, structured cell
> filtering and sorting strategies (with verbatim protocol excerpts), clinical response annotations, and
> direct downloadable matrix URLs for all 352 discovered cohorts across 10 major solid tumor indications.

## 1. Table of Contents & Quick Navigation

- [Executive Summary & Multi-Cancer Statistics](#2-executive-summary--multi-cancer-statistics)
- [Cell Filtering & Technology Distributions](#3-cell-filtering--technology-distributions)
- [Master Comparison Table (All 352 Cohorts)](#4-master-comparison-table-all-352-cohorts)
- [Cohort Directory by Cancer Indication](#5-cohort-directory-by-cancer-indication)
  - [Melanoma Cohorts](#melanoma-cohorts)
  - [NSCLC Cohorts](#nsclc-cohorts)
  - [ccRCC Cohorts](#ccrcc-cohorts)
  - [Bladder Cohorts](#bladder-cohorts)
  - [Breast Cohorts](#breast-cohorts)
  - [CRC Cohorts](#crc-cohorts)
  - [HNSCC Cohorts](#hnscc-cohorts)
  - [Gastric Cohorts](#gastric-cohorts)
  - [HCC Cohorts](#hcc-cohorts)
  - [PDAC Cohorts](#pdac-cohorts)

---

## 2. Executive Summary & Multi-Cancer Statistics

| Indication | Tier | Cohorts | Samples/Patients | Est. Total Cells |
|:---|:---|:---:|:---:|:---:|
| **Bladder** | Tier 2 (Baseline Atlas) | 18 | 720 | 1,463,059 |
| **Breast** | Tier 1 (ICB Response) | 12 | 523 | 1,046,000 |
| **Breast** | Tier 1 (ICB Treated) | 2 | 4 | 8,000 |
| **Breast** | Tier 2 (Baseline Atlas) | 75 | 1,055 | 3,212,955 |
| **CRC** | Tier 1 (ICB Response) | 9 | 286 | 572,000 |
| **CRC** | Tier 1 (ICB Treated) | 10 | 110 | 141,190 |
| **CRC** | Tier 2 (Baseline Atlas) | 37 | 2,746 | 8,014,200 |
| **Gastric** | Tier 1 (ICB Response) | 2 | 145 | 290,000 |
| **Gastric** | Tier 2 (Baseline Atlas) | 14 | 272 | 679,611 |
| **HCC** | Tier 1 (ICB Response) | 12 | 313 | 681,466 |
| **HCC** | Tier 1 (ICB Treated) | 4 | 248 | 496,000 |
| **HCC** | Tier 2 (Baseline Atlas) | 6 | 102 | 204,000 |
| **HNSCC** | Tier 1 (ICB Response) | 5 | 317 | 1,140,399 |
| **HNSCC** | Tier 1 (ICB Treated) | 3 | 75 | 150,000 |
| **HNSCC** | Tier 2 (Baseline Atlas) | 13 | 314 | 1,033,949 |
| **Melanoma** | Tier 1 (ICB Response) | 12 | 280 | 560,000 |
| **Melanoma** | Tier 1 (ICB Treated) | 3 | 186 | 393,941 |
| **Melanoma** | Tier 2 (Baseline Atlas) | 14 | 285 | 903,992 |
| **NSCLC** | Tier 1 (ICB Response) | 13 | 604 | 1,208,000 |
| **NSCLC** | Tier 1 (ICB Treated) | 1 | 6 | 12,000 |
| **NSCLC** | Tier 2 (Baseline Atlas) | 22 | 1,786 | 5,753,221 |
| **PDAC** | Tier 1 (ICB Response) | 6 | 220 | 440,000 |
| **PDAC** | Tier 1 (ICB Treated) | 5 | 136 | 272,000 |
| **PDAC** | Tier 2 (Baseline Atlas) | 10 | 147 | 315,905 |
| **ccRCC** | Tier 1 (ICB Response) | 5 | 56 | 112,000 |
| **ccRCC** | Tier 1 (ICB Treated) | 2 | 17 | 34,000 |
| **ccRCC** | Tier 2 (Baseline Atlas) | 37 | 173 | 1,028,201 |

---

## 3. Cell Filtering & Technology Distributions

### Cell Selection & Filtering Strategies Across Cohorts

| Cell Selection / Isolation Strategy | Cohorts | Est. Total Cells | Description & Scope |
|:---|:---:|:---:|:---|
| **Unselected / Total Single-Cell Suspension** | 286 | 26,648,439 | Unselected single-cell suspension. |
| **FACS-sorted (CD3+/CD8+ T-cell enriched)** | 46 | 2,439,433 | Targeted FACS sorting of T lymphocytes/TILs. |
| **Nuclei Isolation (snRNA-seq)** | 14 | 808,217 | Nuclear lysis from frozen archival tissue (snRNA-seq). |
| **FACS-sorted (CD45+ Immune-enriched)** | 3 | 70,000 | FACS-sorted immune compartment; excludes malignant & stromal cells. |
| **FACS-sorted (EpCAM+ Malignant/Epithelial)** | 1 | 20,000 | Targeted FACS sorting of malignant/epithelial cells. |
| **FACS-sorted (Viability DAPI-/7-AAD- only)** | 1 | 10,000 | Viable gating only (DAPI-/7-AAD-); no lineage bias. |
| **MACS Bead-selected** | 1 | 170,000 | Unselected single-cell suspension. |

### Sequencing Technology & Chemistry Distribution

| Platform / Chemistry | Cohorts | Est. Total Cells |
|:---|:---:|:---:|
| **High-Throughput scRNA-seq** | 187 | 10,418,913 |
| **10x Chromium 3' v3/v3.1** | 69 | 16,964,540 |
| **10x Visium Spatial Transcriptomics** | 53 | 893,859 |
| **Single-Nucleus RNA-seq** | 14 | 660,000 |
| **10x Chromium (3' unspecified)** | 10 | 214,643 |
| **10x Chromium 3' v2** | 9 | 173,194 |
| **10x Chromium 5' (Immune Profiling)** | 7 | 657,836 |
| **Subcellular Spatial Transcriptomics (CosMx/Xenium)** | 2 | 108,000 |
| **BD Rhapsody** | 1 | 75,104 |

---

## 4. Master Comparison Table (All 352 Cohorts)

| # | Accession | Indication | Tier | Technology | Cell Selection / Filtering | Samples | Est. Cells | Response? | Repository | DOI / Reference |
|:---:|:---|:---|:---|:---|:---|:---:|:---:|:---:|:---|:---:|
| 1 | [GSE246613](#gse246613) | **Breast** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 266 | 532,000 | Yes | NCBI GEO | [PMID:38194915](https://doi.org/38194915) |
| 2 | [GSE274141](#gse274141) | **Breast** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 54 | 108,000 | Yes | NCBI GEO | [PMID:39955556](https://doi.org/39955556) |
| 3 | [GSE222859](#gse222859) | **Breast** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 51 | 102,000 | Yes | NCBI GEO | [PMID:37248301](https://doi.org/37248301) |
| 4 | [GSE274139](#gse274139) | **Breast** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 35 | 70,000 | Yes | NCBI GEO | [PMID:39955556](https://doi.org/39955556) |
| 5 | [GSE300475](#gse300475) | **Breast** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 32 | 64,000 | Yes | NCBI GEO | [PMID:40610460](https://doi.org/40610460) |
| 6 | [GSE212707](#gse212707) | **Breast** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 21 | 42,000 | Yes | NCBI GEO | GEO / CZI |
| 7 | [GSE262288](#gse262288) | **Breast** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 19 | 38,000 | Yes | NCBI GEO | [PMID:39955556](https://doi.org/39955556) |
| 8 | [GSE303346](#gse303346) | **Breast** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 16 | 32,000 | Yes | NCBI GEO | GEO / CZI |
| 9 | [GSE254991](#gse254991) | **Breast** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 12 | 24,000 | Yes | NCBI GEO | [PMID:40064918](https://doi.org/40064918) |
| 10 | [GSE199219](#gse199219) | **Breast** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 10 | 20,000 | Yes | NCBI GEO | GEO / CZI |
| 11 | [GSE332708](#gse332708) | **Breast** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 4 | 8,000 | Yes | NCBI GEO | GEO / CZI |
| 12 | [GSE331487](#gse331487) | **Breast** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 3 | 6,000 | Yes | NCBI GEO | GEO / CZI |
| 13 | [GSE236581](#gse236581) | **CRC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 169 | 338,000 | Yes | NCBI GEO | [PMID:38981439](https://doi.org/38981439) |
| 14 | [GSE205506](#gse205506) | **CRC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 40 | 80,000 | Yes | NCBI GEO | [PMID:37172580](https://doi.org/37172580) |
| 15 | [GSE299651](#gse299651) | **CRC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 20 | 40,000 | Yes | NCBI GEO | [PMID:41555096](https://doi.org/41555096) |
| 16 | [GSE146771](#gse146771) | **CRC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 20 | 40,000 | Yes | NCBI GEO | [PMID:32302573](https://doi.org/32302573) |
| 17 | [GSE164522](#gse164522) | **CRC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD45+ Immune-enriched) | 17 | 34,000 | Yes | NCBI GEO | [PMID:35303421](https://doi.org/35303421) |
| 18 | [GSE235917](#gse235917) | **CRC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 8 | 16,000 | Yes | NCBI GEO | [PMID:37563240](https://doi.org/37563240) |
| 19 | [GSE278406](#gse278406) | **CRC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 6 | 12,000 | Yes | NCBI GEO | [PMID:39793571](https://doi.org/39793571) |
| 20 | [GSE274321](#gse274321) | **CRC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 3 | 6,000 | Yes | NCBI GEO | [PMID:41484106](https://doi.org/41484106) |
| 21 | [GSE309346](#gse309346) | **CRC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 3 | 6,000 | Yes | NCBI GEO | [PMID:41817574](https://doi.org/41817574) |
| 22 | [GSE270680](#gse270680) | **Gastric** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 77 | 154,000 | Yes | NCBI GEO | [PMID:41593079](https://doi.org/41593079) |
| 23 | [GSE239676](#gse239676) | **Gastric** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 68 | 136,000 | Yes | NCBI GEO | [PMID:39097198](https://doi.org/39097198) |
| 24 | [GSE313642](#gse313642) | **HCC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 194 | 388,000 | Yes | NCBI GEO | [PMID:41831609](https://doi.org/41831609) |
| 25 | [GSE245906](#gse245906) | **HCC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 20 | 40,000 | Yes | NCBI GEO | [PMID:38350444](https://doi.org/38350444) |
| 26 | [GSE319709](#gse319709) | **HCC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 18 | 36,000 | Yes | NCBI GEO | GEO / CZI |
| 27 | [GSE272347](#gse272347) | **HCC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 17 | 34,000 | Yes | NCBI GEO | [PMID:39455893](https://doi.org/39455893) |
| 28 | [GSE318418](#gse318418) | **HCC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 16 | 32,000 | Yes | NCBI GEO | [PMID:41831609](https://doi.org/41831609) |
| 29 | [GSE299340](#gse299340) | **HCC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 10 | 20,000 | Yes | NCBI GEO | [PMID:40721811](https://doi.org/40721811) |
| 30 | [GSE282343](#gse282343) | **HCC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 8 | 71,466 | Yes | NCBI GEO | [PMID:40634523](https://doi.org/40634523) |
| 31 | [GSE255830](#gse255830) | **HCC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 8 | 16,000 | Yes | NCBI GEO | [PMID:38584166](https://doi.org/38584166) |
| 32 | [GSE281110](#gse281110) | **HCC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 7 | 14,000 | Yes | NCBI GEO | [PMID:39934824](https://doi.org/39934824) |
| 33 | [GSE233405](#gse233405) | **HCC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 6 | 12,000 | Yes | NCBI GEO | [PMID:37342338](https://doi.org/37342338) |
| 34 | [GSE278324](#gse278324) | **HCC** | Tier 1 (ICB Response) | Single-Nucleus RNA-seq | Nuclei Isolation (snRNA-seq) | 6 | 12,000 | Yes | NCBI GEO | GEO / CZI |
| 35 | [GSE265770](#gse265770) | **HCC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 3 | 6,000 | Yes | NCBI GEO | [PMID:39318093](https://doi.org/39318093) |
| 36 | [GSE200996](#gse200996) | **HNSCC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 204 | 408,000 | Yes | NCBI GEO | [PMID:35803260](https://doi.org/35803260) |
| 37 | [GSE301741](#gse301741) | **HNSCC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 58 | 116,000 | Yes | NCBI GEO | [PMID:41923630](https://doi.org/41923630) |
| 38 | [GSE296954](#gse296954) | **HNSCC** | Tier 1 (ICB Response) | 10x Chromium 5' (Immune Profiling) | FACS-sorted (CD3+/CD8+ T-cell enriched) | 44 | 88,000 | Yes | NCBI GEO | GEO / CZI |
| 39 | [GSE301720](#gse301720) | **HNSCC** | Tier 1 (ICB Response) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 8 | 522,399 | Yes | NCBI GEO | [PMID:41837744](https://doi.org/41837744) |
| 40 | [GSE247582](#gse247582) | **HNSCC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 3 | 6,000 | Yes | NCBI GEO | [PMID:38077321](https://doi.org/38077321) |
| 41 | [GSE344166](#gse344166) | **Melanoma** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 116 | 232,000 | Yes | NCBI GEO | GEO / CZI |
| 42 | [GSE286410](#gse286410) | **Melanoma** | Tier 1 (ICB Response) | Single-Nucleus RNA-seq | Unselected / Total Single-Cell Suspension | 41 | 82,000 | Yes | NCBI GEO | [PMID:41168392](https://doi.org/41168392) |
| 43 | [GSE218429](#gse218429) | **Melanoma** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 35 | 70,000 | Yes | NCBI GEO | [PMID:36864091](https://doi.org/36864091) |
| 44 | [GSE198265](#gse198265) | **Melanoma** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 22 | 44,000 | Yes | NCBI GEO | [PMID:35413271](https://doi.org/35413271) |
| 45 | [GSE242477](#gse242477) | **Melanoma** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 17 | 34,000 | Yes | NCBI GEO | [PMID:39630887](https://doi.org/39630887) |
| 46 | [GSE303948](#gse303948) | **Melanoma** | Tier 1 (ICB Response) | Single-Nucleus RNA-seq | Unselected / Total Single-Cell Suspension | 14 | 28,000 | Yes | NCBI GEO | [PMID:40890106](https://doi.org/40890106) |
| 47 | [GSE294273](#gse294273) | **Melanoma** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 12 | 24,000 | Yes | NCBI GEO | [PMID:42277002](https://doi.org/42277002) |
| 48 | [GSE320040](#gse320040) | **Melanoma** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 7 | 14,000 | Yes | NCBI GEO | [PMID:42092150](https://doi.org/42092150) |
| 49 | [GSE210963](#gse210963) | **Melanoma** | Tier 1 (ICB Response) | 10x Chromium (3' unspecified) | Unselected / Total Single-Cell Suspension | 6 | 12,000 | Yes | NCBI GEO | [PMID:38015751](https://doi.org/38015751) |
| 50 | [GSE256291](#gse256291) | **Melanoma** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 4 | 8,000 | Yes | NCBI GEO | [PMID:39266214](https://doi.org/39266214) |
| 51 | [GSE244983](#gse244983) | **Melanoma** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 4 | 8,000 | Yes | NCBI GEO | [PMID:38594286](https://doi.org/38594286) |
| 52 | [GSE300446](#gse300446) | **Melanoma** | Tier 1 (ICB Response) | Single-Nucleus RNA-seq | Nuclei Isolation (snRNA-seq) | 2 | 4,000 | Yes | NCBI GEO | GEO / CZI |
| 53 | [GSE243013](#gse243013) | **NSCLC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 243 | 486,000 | Yes | NCBI GEO | [PMID:40147443](https://doi.org/40147443) |
| 54 | [GSE241934](#gse241934) | **NSCLC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 88 | 176,000 | Yes | NCBI GEO | [PMID:38897205](https://doi.org/38897205) |

*Full 352-cohort details available in docs/comprehensive_solid_tumor_sc_datasets_catalog.md*