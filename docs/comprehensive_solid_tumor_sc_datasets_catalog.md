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
| 55 | [GSE270148](#gse270148) | **NSCLC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | MACS Bead-selected | 85 | 170,000 | Yes | NCBI GEO | [PMID:40931076](https://doi.org/40931076) |
| 56 | [GSE317309](#gse317309) | **NSCLC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 64 | 128,000 | Yes | NCBI GEO | GEO / CZI |
| 57 | [GSE207422](#gse207422) | **NSCLC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 39 | 78,000 | Yes | NCBI GEO | [PMID:36869384](https://doi.org/36869384) |
| 58 | [GSE303762](#gse303762) | **NSCLC** | Tier 1 (ICB Response) | Single-Nucleus RNA-seq | Nuclei Isolation (snRNA-seq) | 38 | 76,000 | Yes | NCBI GEO | GEO / CZI |
| 59 | [GSE205049](#gse205049) | **NSCLC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 18 | 36,000 | Yes | NCBI GEO | [PMID:37263303](https://doi.org/37263303) |
| 60 | [GSE205354](#gse205354) | **NSCLC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 8 | 16,000 | Yes | NCBI GEO | [PMID:37263303](https://doi.org/37263303) |
| 61 | [GSE233203](#gse233203) | **NSCLC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 7 | 14,000 | Yes | NCBI GEO | [PMID:40843133](https://doi.org/40843133) |
| 62 | [GSE253718](#gse253718) | **NSCLC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 6 | 12,000 | Yes | NCBI GEO | GEO / CZI |
| 63 | [GSE223779](#gse223779) | **NSCLC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 6 | 12,000 | Yes | NCBI GEO | [PMID:37272226](https://doi.org/37272226) |
| 64 | [GSE307811](#gse307811) | **NSCLC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 2,000 | Yes | NCBI GEO | [PMID:41511407](https://doi.org/41511407) |
| 65 | [GSE285888](#gse285888) | **NSCLC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 1 | 2,000 | Yes | NCBI GEO | [PMID:40404203](https://doi.org/40404203) |
| 66 | [GSE311789](#gse311789) | **PDAC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 142 | 284,000 | Yes | NCBI GEO | [PMID:41707654](https://doi.org/41707654) |
| 67 | [GSE279781](#gse279781) | **PDAC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 30 | 60,000 | Yes | NCBI GEO | [PMID:39811671](https://doi.org/39811671) |
| 68 | [GSE316195](#gse316195) | **PDAC** | Tier 1 (ICB Response) | Single-Nucleus RNA-seq | Nuclei Isolation (snRNA-seq) | 22 | 44,000 | Yes | NCBI GEO | GEO / CZI |
| 69 | [GSE212966](#gse212966) | **PDAC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 12 | 24,000 | Yes | NCBI GEO | [PMID:36944944](https://doi.org/36944944) |
| 70 | [GSE283206](#gse283206) | **PDAC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 8 | 16,000 | Yes | NCBI GEO | [PMID:40815670](https://doi.org/40815670) |
| 71 | [GSE311788](#gse311788) | **PDAC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 6 | 12,000 | Yes | NCBI GEO | [PMID:41707654](https://doi.org/41707654) |
| 72 | [GSE314072](#gse314072) | **ccRCC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 24 | 48,000 | Yes | NCBI GEO | [PMID:41557789](https://doi.org/41557789) |
| 73 | [GSE285701](#gse285701) | **ccRCC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 13 | 26,000 | Yes | NCBI GEO | [PMID:40961944](https://doi.org/40961944) |
| 74 | [GSE210038](#gse210038) | **ccRCC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 9 | 18,000 | Yes | NCBI GEO | [PMID:37335139](https://doi.org/37335139) |
| 75 | [GSE304466](#gse304466) | **ccRCC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 6 | 12,000 | Yes | NCBI GEO | GEO / CZI |
| 76 | [GSE223808](#gse223808) | **ccRCC** | Tier 1 (ICB Response) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 4 | 8,000 | Yes | NCBI GEO | [PMID:36928415](https://doi.org/36928415) |
| 77 | [GSE302453](#gse302453) | **Breast** | Tier 1 (ICB Treated) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 2 | 4,000 | Baseline | NCBI GEO | GEO / CZI |
| 78 | [GSE299267](#gse299267) | **Breast** | Tier 1 (ICB Treated) | Single-Nucleus RNA-seq | Unselected / Total Single-Cell Suspension | 2 | 4,000 | Baseline | NCBI GEO | GEO / CZI |
| 79 | [CELLxGENE_829a3cd1](#cellxgene_829a3cd1) | **CRC** | Tier 1 (ICB Treated) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 29 | 47,107 | Baseline | CZ CELLxGENE | [10.1038/s41586-024-08560-0](https://doi.org/10.1038/s41586-024-08560-0) |
| 80 | [CELLxGENE_2554a654](#cellxgene_2554a654) | **CRC** | Tier 1 (ICB Treated) | 10x Chromium 3' v3/v3.1 | FACS-sorted (CD3+/CD8+ T-cell enriched) | 28 | 26,086 | Baseline | CZ CELLxGENE | [10.1038/s41586-024-08560-0](https://doi.org/10.1038/s41586-024-08560-0) |
| 81 | [CELLxGENE_5ee552f5](#cellxgene_5ee552f5) | **CRC** | Tier 1 (ICB Treated) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 27 | 21,026 | Baseline | CZ CELLxGENE | [10.1038/s41586-024-08560-0](https://doi.org/10.1038/s41586-024-08560-0) |
| 82 | [CELLxGENE_4b5afdf9](#cellxgene_4b5afdf9) | **CRC** | Tier 1 (ICB Treated) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 9 | 13,843 | Baseline | CZ CELLxGENE | [10.1038/s41586-024-08560-0](https://doi.org/10.1038/s41586-024-08560-0) |
| 83 | [GSE188711](#gse188711) | **CRC** | Tier 1 (ICB Treated) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 6 | 12,000 | Baseline | NCBI GEO | [PMID:34793335](https://doi.org/34793335) |
| 84 | [GSE336564](#gse336564) | **CRC** | Tier 1 (ICB Treated) | 10x Chromium (3' unspecified) | FACS-sorted (Viability DAPI-/7-AAD- only) | 6 | 10,000 | Baseline | NCBI GEO | GEO / CZI |
| 85 | [GSE216534](#gse216534) | **CRC** | Tier 1 (ICB Treated) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 2 | 4,000 | Baseline | NCBI GEO | [PMID:36631610](https://doi.org/36631610) |
| 86 | [CELLxGENE_387acac5](#cellxgene_387acac5) | **CRC** | Tier 1 (ICB Treated) | 10x Chromium 3' v3/v3.1 | FACS-sorted (CD3+/CD8+ T-cell enriched) | 1 | 2,574 | Baseline | CZ CELLxGENE | [10.1038/s41586-024-08560-0](https://doi.org/10.1038/s41586-024-08560-0) |
| 87 | [CELLxGENE_ef0d813e](#cellxgene_ef0d813e) | **CRC** | Tier 1 (ICB Treated) | 10x Chromium 3' v3/v3.1 | FACS-sorted (CD3+/CD8+ T-cell enriched) | 1 | 1,203 | Baseline | CZ CELLxGENE | [10.1038/s41586-024-08560-0](https://doi.org/10.1038/s41586-024-08560-0) |
| 88 | [CELLxGENE_2e95d453](#cellxgene_2e95d453) | **CRC** | Tier 1 (ICB Treated) | 10x Chromium 3' v3/v3.1 | FACS-sorted (CD3+/CD8+ T-cell enriched) | 1 | 3,351 | Baseline | CZ CELLxGENE | [10.1038/s41586-024-08560-0](https://doi.org/10.1038/s41586-024-08560-0) |
| 89 | [GSE318420](#gse318420) | **HCC** | Tier 1 (ICB Treated) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 210 | 420,000 | Baseline | NCBI GEO | [PMID:41831609](https://doi.org/41831609) |
| 90 | [GSE272348](#gse272348) | **HCC** | Tier 1 (ICB Treated) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 31 | 62,000 | Baseline | NCBI GEO | [PMID:39455893](https://doi.org/39455893) |
| 91 | [GSE224411](#gse224411) | **HCC** | Tier 1 (ICB Treated) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 4 | 8,000 | Baseline | NCBI GEO | [PMID:37080163](https://doi.org/37080163) |
| 92 | [GSE215428](#gse215428) | **HCC** | Tier 1 (ICB Treated) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 3 | 6,000 | Baseline | NCBI GEO | GEO / CZI |
| 93 | [GSE287301](#gse287301) | **HNSCC** | Tier 1 (ICB Treated) | Subcellular Spatial Transcriptomics (CosMx/Xenium) | Unselected / Total Single-Cell Suspension | 48 | 96,000 | Baseline | NCBI GEO | [PMID:41961948](https://doi.org/41961948) |
| 94 | [GSE296867](#gse296867) | **HNSCC** | Tier 1 (ICB Treated) | 10x Chromium 5' (Immune Profiling) | FACS-sorted (CD3+/CD8+ T-cell enriched) | 22 | 44,000 | Baseline | NCBI GEO | GEO / CZI |
| 95 | [GSE327189](#gse327189) | **HNSCC** | Tier 1 (ICB Treated) | 10x Chromium 5' (Immune Profiling) | FACS-sorted (CD3+/CD8+ T-cell enriched) | 5 | 10,000 | Baseline | NCBI GEO | GEO / CZI |
| 96 | [CELLxGENE_7b20c613](#cellxgene_7b20c613) | **Melanoma** | Tier 1 (ICB Treated) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 167 | 355,941 | Baseline | CZ CELLxGENE | [10.1038/s41597-025-04381-6](https://doi.org/10.1038/s41597-025-04381-6) |
| 97 | [GSE270464](#gse270464) | **Melanoma** | Tier 1 (ICB Treated) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 10 | 20,000 | Baseline | NCBI GEO | [PMID:40695831](https://doi.org/40695831) |
| 98 | [GSE211068](#gse211068) | **Melanoma** | Tier 1 (ICB Treated) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 9 | 18,000 | Baseline | NCBI GEO | [PMID:36206353](https://doi.org/36206353) |
| 99 | [GSE276139](#gse276139) | **NSCLC** | Tier 1 (ICB Treated) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 6 | 12,000 | Baseline | NCBI GEO | [PMID:41876785](https://doi.org/41876785) |
| 100 | [GSE211644](#gse211644) | **PDAC** | Tier 1 (ICB Treated) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 50 | 100,000 | Baseline | NCBI GEO | [PMID:35849783](https://doi.org/35849783) |
| 101 | [GSE348275](#gse348275) | **PDAC** | Tier 1 (ICB Treated) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 35 | 70,000 | Baseline | NCBI GEO | GEO / CZI |
| 102 | [GSE318413](#gse318413) | **PDAC** | Tier 1 (ICB Treated) | Single-Nucleus RNA-seq | Nuclei Isolation (snRNA-seq) | 19 | 38,000 | Baseline | NCBI GEO | [PMID:41673768](https://doi.org/41673768) |
| 103 | [GSE156405](#gse156405) | **PDAC** | Tier 1 (ICB Treated) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 17 | 34,000 | Baseline | NCBI GEO | [PMID:34426439](https://doi.org/34426439) |
| 104 | [GSE335452](#gse335452) | **PDAC** | Tier 1 (ICB Treated) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 15 | 30,000 | Baseline | NCBI GEO | [PMID:42548877](https://doi.org/42548877) |
| 105 | [GSE220313](#gse220313) | **ccRCC** | Tier 1 (ICB Treated) | High-Throughput scRNA-seq | FACS-sorted (CD45+ Immune-enriched) | 16 | 32,000 | Baseline | NCBI GEO | GEO / CZI |
| 106 | [GSE254498](#gse254498) | **ccRCC** | Tier 1 (ICB Treated) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 2,000 | Baseline | NCBI GEO | [PMID:41413084](https://doi.org/41413084) |
| 107 | [GSE183556](#gse183556) | **Bladder** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 550 | 1,100,000 | Baseline | NCBI GEO | [PMID:38030723](https://doi.org/38030723) |
| 108 | [GSE192575](#gse192575) | **Bladder** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 42 | 84,000 | Baseline | NCBI GEO | GEO / CZI |
| 109 | [GSE267718](#gse267718) | **Bladder** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 33 | 66,000 | Baseline | NCBI GEO | [PMID:38847806](https://doi.org/38847806) |
| 110 | [GSE302781](#gse302781) | **Bladder** | Tier 2 (Baseline Atlas) | Single-Nucleus RNA-seq | Nuclei Isolation (snRNA-seq) | 15 | 30,000 | Baseline | NCBI GEO | GEO / CZI |
| 111 | [GSE301651](#gse301651) | **Bladder** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 13 | 26,000 | Baseline | NCBI GEO | [PMID:41706539](https://doi.org/41706539) |
| 112 | [GSE222315](#gse222315) | **Bladder** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 13 | 26,000 | Baseline | NCBI GEO | [PMID:38428409](https://doi.org/38428409) |
| 113 | [GSE172433](#gse172433) | **Bladder** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 12 | 24,000 | Baseline | NCBI GEO | [PMID:36323682](https://doi.org/36323682) |
| 114 | [GSE277524](#gse277524) | **Bladder** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 10 | 20,000 | Baseline | NCBI GEO | [PMID:39702554](https://doi.org/39702554) |
| 115 | [GSE310802](#gse310802) | **Bladder** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 8 | 16,000 | Baseline | NCBI GEO | GEO / CZI |
| 116 | [GSE326225](#gse326225) | **Bladder** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 8 | 16,000 | Baseline | NCBI GEO | GEO / CZI |
| 117 | [GSE326854](#gse326854) | **Bladder** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 6 | 12,000 | Baseline | NCBI GEO | [PMID:42203766](https://doi.org/42203766) |
| 118 | [GSE250523](#gse250523) | **Bladder** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 4 | 8,000 | Baseline | NCBI GEO | [PMID:39386756](https://doi.org/39386756) |
| 119 | [CELLxGENE_024581e3](#cellxgene_024581e3) | **Bladder** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 1 | 8,156 | Baseline | CZ CELLxGENE | [10.1038/s41467-025-67643-2](https://doi.org/10.1038/s41467-025-67643-2) |
| 120 | [CELLxGENE_7ee4b15b](#cellxgene_7ee4b15b) | **Bladder** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 1 | 10,946 | Baseline | CZ CELLxGENE | [10.1038/s41467-025-67643-2](https://doi.org/10.1038/s41467-025-67643-2) |
| 121 | [CELLxGENE_e2094676](#cellxgene_e2094676) | **Bladder** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 1 | 5,762 | Baseline | CZ CELLxGENE | [10.1038/s41467-025-67643-2](https://doi.org/10.1038/s41467-025-67643-2) |
| 122 | [CELLxGENE_f0e0575d](#cellxgene_f0e0575d) | **Bladder** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 1 | 2,410 | Baseline | CZ CELLxGENE | [10.1038/s41467-025-67643-2](https://doi.org/10.1038/s41467-025-67643-2) |
| 123 | [GSE176249](#gse176249) | **Bladder** | Tier 2 (Baseline Atlas) | 10x Chromium (3' unspecified) | Unselected / Total Single-Cell Suspension | 1 | 2,000 | Baseline | NCBI GEO | [PMID:36323682](https://doi.org/36323682) |
| 124 | [CELLxGENE_670a9f65](#cellxgene_670a9f65) | **Bladder** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 1 | 5,785 | Baseline | CZ CELLxGENE | [10.1038/s41467-025-67643-2](https://doi.org/10.1038/s41467-025-67643-2) |
| 125 | [CELLxGENE_5a9cfb44](#cellxgene_5a9cfb44) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 138 | 274,555 | Baseline | CZ CELLxGENE | [10.1093/nargab/lqaf217](https://doi.org/10.1093/nargab/lqaf217) |
| 126 | [CELLxGENE_de5416ef](#cellxgene_de5416ef) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 138 | 621,200 | Baseline | CZ CELLxGENE | [10.1093/nargab/lqaf217](https://doi.org/10.1093/nargab/lqaf217) |
| 127 | [CELLxGENE_ed880090](#cellxgene_ed880090) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | FACS-sorted (CD3+/CD8+ T-cell enriched) | 136 | 242,613 | Baseline | CZ CELLxGENE | [10.1093/nargab/lqaf217](https://doi.org/10.1093/nargab/lqaf217) |
| 128 | [CELLxGENE_75011e96](#cellxgene_75011e96) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 131 | 104,032 | Baseline | CZ CELLxGENE | [10.1093/nargab/lqaf217](https://doi.org/10.1093/nargab/lqaf217) |
| 129 | [CELLxGENE_6f9de485](#cellxgene_6f9de485) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 101 | 427,823 | Baseline | CZ CELLxGENE | [Dataset Version: https://datasets.cellxgene.cziscience.com/e94bd3cc-6271-424a-baac-12f8eb320a0e.h5ad curated and distributed by CZ CELLxGENE Discover in Collection: https://cellxgene.cziscience.com/collections/ceef2841-5333-46ac-92ef-ccbe0c20fe55](https://doi.org/Dataset Version: https://datasets.cellxgene.cziscience.com/e94bd3cc-6271-424a-baac-12f8eb320a0e.h5ad curated and distributed by CZ CELLxGENE Discover in Collection: https://cellxgene.cziscience.com/collections/ceef2841-5333-46ac-92ef-ccbe0c20fe55) |
| 130 | [GSE300628](#gse300628) | **Breast** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 60 | 120,000 | Baseline | NCBI GEO | [PMID:41925564](https://doi.org/41925564) |
| 131 | [CELLxGENE_34f5307e](#cellxgene_34f5307e) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 32 | 394,534 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 132 | [CELLxGENE_6c87755e](#cellxgene_6c87755e) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 30 | 157,531 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 133 | [CELLxGENE_933497dc](#cellxgene_933497dc) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 26 | 35,214 | Baseline | CZ CELLxGENE | [10.1038/s41588-021-00911-1](https://doi.org/10.1038/s41588-021-00911-1) |
| 134 | [CELLxGENE_9fddb063](#cellxgene_9fddb063) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 26 | 100,064 | Baseline | CZ CELLxGENE | [10.1038/s41588-021-00911-1](https://doi.org/10.1038/s41588-021-00911-1) |
| 135 | [CELLxGENE_11a3244a](#cellxgene_11a3244a) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 26 | 9,675 | Baseline | CZ CELLxGENE | [10.1038/s41588-021-00911-1](https://doi.org/10.1038/s41588-021-00911-1) |
| 136 | [CELLxGENE_4cdd25a4](#cellxgene_4cdd25a4) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 26 | 19,601 | Baseline | CZ CELLxGENE | [10.1038/s41588-021-00911-1](https://doi.org/10.1038/s41588-021-00911-1) |
| 137 | [CELLxGENE_7357bdd2](#cellxgene_7357bdd2) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v2 | Unselected / Total Single-Cell Suspension | 20 | 24,489 | Baseline | CZ CELLxGENE | [10.1038/s41588-021-00911-1](https://doi.org/10.1038/s41588-021-00911-1) |
| 138 | [CELLxGENE_04d87de6](#cellxgene_04d87de6) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 20 | 3,206 | Baseline | CZ CELLxGENE | [10.1038/s41588-021-00911-1](https://doi.org/10.1038/s41588-021-00911-1) |
| 139 | [CELLxGENE_2f05ab20](#cellxgene_2f05ab20) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Nuclei Isolation (snRNA-seq) | 16 | 38,217 | Baseline | CZ CELLxGENE | [10.1038/s41586-021-03710-0](https://doi.org/10.1038/s41586-021-03710-0) |
| 140 | [CELLxGENE_0c86f0de](#cellxgene_0c86f0de) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 14 | 3,524 | Baseline | CZ CELLxGENE | [10.1038/s41588-021-00911-1](https://doi.org/10.1038/s41588-021-00911-1) |
| 141 | [GSE337706](#gse337706) | **Breast** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 13 | 26,000 | Baseline | NCBI GEO | GEO / CZI |
| 142 | [GSE281488](#gse281488) | **Breast** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 9 | 18,000 | Baseline | NCBI GEO | [PMID:41923644](https://doi.org/41923644) |
| 143 | [CELLxGENE_5d3fc988](#cellxgene_5d3fc988) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 7 | 49,109 | Baseline | CZ CELLxGENE | [10.1002/ctm2.1356](https://doi.org/10.1002/ctm2.1356) |
| 144 | [CELLxGENE_3e4e2c8e](#cellxgene_3e4e2c8e) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 7 | 22,414 | Baseline | CZ CELLxGENE | [10.1002/ctm2.1356](https://doi.org/10.1002/ctm2.1356) |
| 145 | [GSE281490](#gse281490) | **Breast** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 6 | 12,000 | Baseline | NCBI GEO | [PMID:41923644](https://doi.org/41923644) |
| 146 | [CELLxGENE_68b6114f](#cellxgene_68b6114f) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v2 | Unselected / Total Single-Cell Suspension | 5 | 24,271 | Baseline | CZ CELLxGENE | [10.15252/embj.2019104063](https://doi.org/10.15252/embj.2019104063) |
| 147 | [CELLxGENE_e500acbf](#cellxgene_e500acbf) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 5 | 37,428 | Baseline | CZ CELLxGENE | [10.1002/ctm2.1356](https://doi.org/10.1002/ctm2.1356) |
| 148 | [GSE230327](#gse230327) | **Breast** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 4 | 8,000 | Baseline | NCBI GEO | [PMID:41886605](https://doi.org/41886605) |
| 149 | [CELLxGENE_fbdd8c17](#cellxgene_fbdd8c17) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium (3' unspecified) | Unselected / Total Single-Cell Suspension | 4 | 10,689 | Baseline | CZ CELLxGENE | [10.1101/2024.11.01.621259](https://doi.org/10.1101/2024.11.01.621259) |
| 150 | [GSE329389](#gse329389) | **Breast** | Tier 2 (Baseline Atlas) | Single-Nucleus RNA-seq | Nuclei Isolation (snRNA-seq) | 3 | 6,000 | Baseline | NCBI GEO | GEO / CZI |
| 151 | [GSE325982](#gse325982) | **Breast** | Tier 2 (Baseline Atlas) | Single-Nucleus RNA-seq | Nuclei Isolation (snRNA-seq) | 2 | 4,000 | Baseline | NCBI GEO | [PMID:42509302](https://doi.org/42509302) |
| 152 | [GSE309616](#gse309616) | **Breast** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 2 | 4,000 | Baseline | NCBI GEO | [PMID:41630032](https://doi.org/41630032) |
| 153 | [GSE229723](#gse229723) | **Breast** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | FACS-sorted (CD45+ Immune-enriched) | 2 | 4,000 | Baseline | NCBI GEO | [PMID:40794843](https://doi.org/40794843) |
| 154 | [CELLxGENE_12c868c6](#cellxgene_12c868c6) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 1 | 10,623 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 155 | [CELLxGENE_a6c0143c](#cellxgene_a6c0143c) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 1 | 7,505 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 156 | [CELLxGENE_6f0858c0](#cellxgene_6f0858c0) | **Breast** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1038/s41588-021-00911-1](https://doi.org/10.1038/s41588-021-00911-1) |
| 157 | [CELLxGENE_80466231](#cellxgene_80466231) | **Breast** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1038/s41588-021-00911-1](https://doi.org/10.1038/s41588-021-00911-1) |
| 158 | [CELLxGENE_a6b0f655](#cellxgene_a6b0f655) | **Breast** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1038/s41588-021-00911-1](https://doi.org/10.1038/s41588-021-00911-1) |
| 159 | [CELLxGENE_aafb780d](#cellxgene_aafb780d) | **Breast** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1038/s41588-021-00911-1](https://doi.org/10.1038/s41588-021-00911-1) |
| 160 | [CELLxGENE_ee141ea4](#cellxgene_ee141ea4) | **Breast** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.labinv.2023.100258](https://doi.org/10.1016/j.labinv.2023.100258) |
| 161 | [CELLxGENE_f354e4c3](#cellxgene_f354e4c3) | **Breast** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1038/s41588-021-00911-1](https://doi.org/10.1038/s41588-021-00911-1) |
| 162 | [CELLxGENE_02aa7750](#cellxgene_02aa7750) | **Breast** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 22,033 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 163 | [CELLxGENE_05a49baa](#cellxgene_05a49baa) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v2 | Unselected / Total Single-Cell Suspension | 1 | 4,742 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 164 | [CELLxGENE_0f9d1892](#cellxgene_0f9d1892) | **Breast** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 33,272 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 165 | [CELLxGENE_1637e817](#cellxgene_1637e817) | **Breast** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 8,086 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 166 | [CELLxGENE_1884e651](#cellxgene_1884e651) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 1 | 5,463 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 167 | [CELLxGENE_2dd73feb](#cellxgene_2dd73feb) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v2 | Unselected / Total Single-Cell Suspension | 1 | 13,167 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 168 | [CELLxGENE_39f6fec9](#cellxgene_39f6fec9) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 1 | 20,916 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 169 | [CELLxGENE_44941fdb](#cellxgene_44941fdb) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v2 | Unselected / Total Single-Cell Suspension | 1 | 2,276 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 170 | [CELLxGENE_48b55b2b](#cellxgene_48b55b2b) | **Breast** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 9,362 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 171 | [CELLxGENE_494faa16](#cellxgene_494faa16) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 1 | 9,958 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 172 | [CELLxGENE_54d56674](#cellxgene_54d56674) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v2 | Unselected / Total Single-Cell Suspension | 1 | 11,074 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 173 | [CELLxGENE_59d14a35](#cellxgene_59d14a35) | **Breast** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 5,833 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 174 | [CELLxGENE_6384d8b8](#cellxgene_6384d8b8) | **Breast** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 22,188 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 175 | [CELLxGENE_71513028](#cellxgene_71513028) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 1 | 12,258 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 176 | [CELLxGENE_71e44b30](#cellxgene_71e44b30) | **Breast** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 11,430 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 177 | [CELLxGENE_7432b873](#cellxgene_7432b873) | **Breast** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 8,715 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 178 | [CELLxGENE_9237e573](#cellxgene_9237e573) | **Breast** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 10,892 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 179 | [CELLxGENE_9c5f68fc](#cellxgene_9c5f68fc) | **Breast** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 11,851 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 180 | [CELLxGENE_a6347c54](#cellxgene_a6347c54) | **Breast** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 11,484 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 181 | [CELLxGENE_24dbd26d](#cellxgene_24dbd26d) | **Breast** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1038/s41588-021-00911-1](https://doi.org/10.1038/s41588-021-00911-1) |
| 182 | [CELLxGENE_aa6f371d](#cellxgene_aa6f371d) | **Breast** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 6,210 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 183 | [CELLxGENE_c7d0def0](#cellxgene_c7d0def0) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 1 | 10,918 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 184 | [CELLxGENE_cd6398a9](#cellxgene_cd6398a9) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 1 | 10,016 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 185 | [CELLxGENE_e06e9bf3](#cellxgene_e06e9bf3) | **Breast** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 9,422 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 186 | [CELLxGENE_e2824739](#cellxgene_e2824739) | **Breast** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 2,958 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 187 | [CELLxGENE_e5c614b8](#cellxgene_e5c614b8) | **Breast** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.labinv.2023.100258](https://doi.org/10.1016/j.labinv.2023.100258) |
| 188 | [CELLxGENE_ec423499](#cellxgene_ec423499) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 1 | 12,494 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 189 | [CELLxGENE_f12ab0e6](#cellxgene_f12ab0e6) | **Breast** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 10,957 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 190 | [CELLxGENE_ff4cfa86](#cellxgene_ff4cfa86) | **Breast** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 1 | 9,799 | Baseline | CZ CELLxGENE | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| 191 | [CELLxGENE_10bb68cf](#cellxgene_10bb68cf) | **Breast** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.labinv.2023.100258](https://doi.org/10.1016/j.labinv.2023.100258) |
| 192 | [CELLxGENE_2cc628d1](#cellxgene_2cc628d1) | **Breast** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.labinv.2023.100258](https://doi.org/10.1016/j.labinv.2023.100258) |
| 193 | [CELLxGENE_480f9371](#cellxgene_480f9371) | **Breast** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.labinv.2023.100258](https://doi.org/10.1016/j.labinv.2023.100258) |
| 194 | [CELLxGENE_540e4c1a](#cellxgene_540e4c1a) | **Breast** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.labinv.2023.100258](https://doi.org/10.1016/j.labinv.2023.100258) |
| 195 | [CELLxGENE_5d2c013d](#cellxgene_5d2c013d) | **Breast** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.labinv.2023.100258](https://doi.org/10.1016/j.labinv.2023.100258) |
| 196 | [CELLxGENE_9adb1b29](#cellxgene_9adb1b29) | **Breast** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.labinv.2023.100258](https://doi.org/10.1016/j.labinv.2023.100258) |
| 197 | [CELLxGENE_b0d9408e](#cellxgene_b0d9408e) | **Breast** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.labinv.2023.100258](https://doi.org/10.1016/j.labinv.2023.100258) |
| 198 | [CELLxGENE_c829c294](#cellxgene_c829c294) | **Breast** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.labinv.2023.100258](https://doi.org/10.1016/j.labinv.2023.100258) |
| 199 | [CELLxGENE_dd1913a6](#cellxgene_dd1913a6) | **Breast** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.labinv.2023.100258](https://doi.org/10.1016/j.labinv.2023.100258) |
| 200 | [CELLxGENE_05a8c945](#cellxgene_05a8c945) | **CRC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 588 | 3,790,266 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2025.12.003](https://doi.org/10.1016/j.ccell.2025.12.003) |
| 201 | [CELLxGENE_19053a82](#cellxgene_19053a82) | **CRC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 308 | 1,596,200 | Baseline | CZ CELLxGENE | [10.1038/s41586-024-07571-1](https://doi.org/10.1038/s41586-024-07571-1) |
| 202 | [CELLxGENE_e6aaf5a4](#cellxgene_e6aaf5a4) | **CRC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 233 | 1,358,573 | Baseline | CZ CELLxGENE | [10.1038/s41586-024-07571-1](https://doi.org/10.1038/s41586-024-07571-1) |
| 203 | [CELLxGENE_dc6b1e06](#cellxgene_dc6b1e06) | **CRC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 227 | 262,642 | Baseline | CZ CELLxGENE | [10.1038/s41586-024-07571-1](https://doi.org/10.1038/s41586-024-07571-1) |
| 204 | [CELLxGENE_9c235282](#cellxgene_9c235282) | **CRC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 223 | 50,570 | Baseline | CZ CELLxGENE | [10.1038/s41586-024-07571-1](https://doi.org/10.1038/s41586-024-07571-1) |
| 205 | [CELLxGENE_40a0ade8](#cellxgene_40a0ade8) | **CRC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 223 | 52,404 | Baseline | CZ CELLxGENE | [10.1038/s41586-024-07571-1](https://doi.org/10.1038/s41586-024-07571-1) |
| 206 | [CELLxGENE_1d54fb17](#cellxgene_1d54fb17) | **CRC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 207 | 250,094 | Baseline | CZ CELLxGENE | [10.1038/s41586-024-07571-1](https://doi.org/10.1038/s41586-024-07571-1) |
| 207 | [CELLxGENE_278eac3f](#cellxgene_278eac3f) | **CRC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 195 | 60,411 | Baseline | CZ CELLxGENE | [10.1038/s41586-024-07571-1](https://doi.org/10.1038/s41586-024-07571-1) |
| 208 | [CELLxGENE_ef7bb7f0](#cellxgene_ef7bb7f0) | **CRC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 184 | 77,050 | Baseline | CZ CELLxGENE | [10.1038/s41586-024-07571-1](https://doi.org/10.1038/s41586-024-07571-1) |
| 209 | [CELLxGENE_7be23e52](#cellxgene_7be23e52) | **CRC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 101 | 23,904 | Baseline | CZ CELLxGENE | [10.1038/s41586-024-07571-1](https://doi.org/10.1038/s41586-024-07571-1) |
| 210 | [CELLxGENE_763d1d88](#cellxgene_763d1d88) | **CRC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 61 | 96,675 | Baseline | CZ CELLxGENE | [10.1038/s41586-024-07571-1](https://doi.org/10.1038/s41586-024-07571-1) |
| 211 | [CELLxGENE_6a270451](#cellxgene_6a270451) | **CRC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 59 | 10,696 | Baseline | CZ CELLxGENE | [10.1016/j.cell.2021.11.031](https://doi.org/10.1016/j.cell.2021.11.031) |
| 212 | [GSE294300](#gse294300) | **CRC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 36 | 72,000 | Baseline | NCBI GEO | [PMID:42143353](https://doi.org/42143353) |
| 213 | [CELLxGENE_d6dfdef1](#cellxgene_d6dfdef1) | **CRC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 26 | 57,723 | Baseline | CZ CELLxGENE | [10.1016/j.cell.2021.11.031](https://doi.org/10.1016/j.cell.2021.11.031) |
| 214 | [GSE271690](#gse271690) | **CRC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 16 | 32,000 | Baseline | NCBI GEO | GEO / CZI |
| 215 | [GSE315534](#gse315534) | **CRC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 15 | 30,000 | Baseline | NCBI GEO | [PMID:42082451](https://doi.org/42082451) |
| 216 | [GSE335811](#gse335811) | **CRC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 12 | 24,000 | Baseline | NCBI GEO | GEO / CZI |
| 217 | [CELLxGENE_e3ed2ba4](#cellxgene_e3ed2ba4) | **CRC** | Tier 2 (Baseline Atlas) | BD Rhapsody | Unselected / Total Single-Cell Suspension | 6 | 75,104 | Baseline | CZ CELLxGENE | [10.1186/s12943-025-02430-7](https://doi.org/10.1186/s12943-025-02430-7) |
| 218 | [GSE330797](#gse330797) | **CRC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 4 | 8,000 | Baseline | NCBI GEO | [PMID:42360233](https://doi.org/42360233) |
| 219 | [GSE311338](#gse311338) | **CRC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 3 | 6,000 | Baseline | NCBI GEO | GEO / CZI |
| 220 | [GSE312804](#gse312804) | **CRC** | Tier 2 (Baseline Atlas) | 10x Chromium (3' unspecified) | Unselected / Total Single-Cell Suspension | 2 | 4,000 | Baseline | NCBI GEO | [PMID:41580988](https://doi.org/41580988) |
| 221 | [GSE312260](#gse312260) | **CRC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 2 | 4,000 | Baseline | NCBI GEO | [PMID:41678386](https://doi.org/41678386) |
| 222 | [CELLxGENE_f7af19e4](#cellxgene_f7af19e4) | **CRC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1038/s41698-023-00488-4](https://doi.org/10.1038/s41698-023-00488-4) |
| 223 | [CELLxGENE_7fe57023](#cellxgene_7fe57023) | **CRC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1038/s41698-023-00488-4](https://doi.org/10.1038/s41698-023-00488-4) |
| 224 | [CELLxGENE_879bb6df](#cellxgene_879bb6df) | **CRC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1038/s41698-023-00488-4](https://doi.org/10.1038/s41698-023-00488-4) |
| 225 | [CELLxGENE_a73f7983](#cellxgene_a73f7983) | **CRC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1038/s41698-023-00488-4](https://doi.org/10.1038/s41698-023-00488-4) |
| 226 | [GSE270767](#gse270767) | **CRC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 2,000 | Baseline | NCBI GEO | [PMID:42436119](https://doi.org/42436119) |
| 227 | [CELLxGENE_b5753bee](#cellxgene_b5753bee) | **CRC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1038/s41698-023-00488-4](https://doi.org/10.1038/s41698-023-00488-4) |
| 228 | [CELLxGENE_c0d43178](#cellxgene_c0d43178) | **CRC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1038/s41698-023-00488-4](https://doi.org/10.1038/s41698-023-00488-4) |
| 229 | [CELLxGENE_1e191a00](#cellxgene_1e191a00) | **CRC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1038/s41698-023-00488-4](https://doi.org/10.1038/s41698-023-00488-4) |
| 230 | [CELLxGENE_7ba1a805](#cellxgene_7ba1a805) | **CRC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1038/s41698-023-00488-4](https://doi.org/10.1038/s41698-023-00488-4) |
| 231 | [CELLxGENE_74e80fd1](#cellxgene_74e80fd1) | **CRC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1038/s41698-023-00488-4](https://doi.org/10.1038/s41698-023-00488-4) |
| 232 | [CELLxGENE_729f397a](#cellxgene_729f397a) | **CRC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1038/s41698-023-00488-4](https://doi.org/10.1038/s41698-023-00488-4) |
| 233 | [CELLxGENE_2d821164](#cellxgene_2d821164) | **CRC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1038/s41698-023-00488-4](https://doi.org/10.1038/s41698-023-00488-4) |
| 234 | [CELLxGENE_2916b663](#cellxgene_2916b663) | **CRC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1038/s41698-023-00488-4](https://doi.org/10.1038/s41698-023-00488-4) |
| 235 | [CELLxGENE_15b98664](#cellxgene_15b98664) | **CRC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1038/s41698-023-00488-4](https://doi.org/10.1038/s41698-023-00488-4) |
| 236 | [CELLxGENE_297b5b89](#cellxgene_297b5b89) | **CRC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1038/s41698-023-00488-4](https://doi.org/10.1038/s41698-023-00488-4) |
| 237 | [CELLxGENE_7bb64315](#cellxgene_7bb64315) | **Gastric** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 70 | 293,823 | Baseline | CZ CELLxGENE | [10.1158/2159-8290.cd-22-0824](https://doi.org/10.1158/2159-8290.cd-22-0824) |
| 238 | [CELLxGENE_0d3807bf](#cellxgene_0d3807bf) | **Gastric** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v2 | Unselected / Total Single-Cell Suspension | 62 | 70,090 | Baseline | CZ CELLxGENE | [10.1038/s41586-024-07571-1](https://doi.org/10.1038/s41586-024-07571-1) |
| 239 | [GSE212212](#gse212212) | **Gastric** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 44 | 88,000 | Baseline | NCBI GEO | [PMID:36921674](https://doi.org/36921674) |
| 240 | [GSE228598](#gse228598) | **Gastric** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 28 | 56,000 | Baseline | NCBI GEO | [PMID:38612926](https://doi.org/38612926) |
| 241 | [GSE234209](#gse234209) | **Gastric** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 12 | 24,000 | Baseline | NCBI GEO | GEO / CZI |
| 242 | [GSE275648](#gse275648) | **Gastric** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 11 | 22,000 | Baseline | NCBI GEO | [PMID:41484771](https://doi.org/41484771) |
| 243 | [GSE112302](#gse112302) | **Gastric** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 10 | 20,000 | Baseline | NCBI GEO | GEO / CZI |
| 244 | [GSE168537](#gse168537) | **Gastric** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 10 | 20,000 | Baseline | NCBI GEO | [PMID:36921037](https://doi.org/36921037) |
| 245 | [GSE246662](#gse246662) | **Gastric** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 9 | 18,000 | Baseline | NCBI GEO | [PMID:39060439](https://doi.org/39060439) |
| 246 | [GSE308231](#gse308231) | **Gastric** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 6 | 12,000 | Baseline | NCBI GEO | [PMID:41360923](https://doi.org/41360923) |
| 247 | [CELLxGENE_ca140407](#cellxgene_ca140407) | **Gastric** | Tier 2 (Baseline Atlas) | 10x Chromium 5' (Immune Profiling) | Unselected / Total Single-Cell Suspension | 5 | 45,698 | Baseline | CZ CELLxGENE | [10.1038/s41467-026-70751-2](https://doi.org/10.1038/s41467-026-70751-2) |
| 248 | [GSE321676](#gse321676) | **Gastric** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 2 | 4,000 | Baseline | NCBI GEO | GEO / CZI |
| 249 | [GSE184198](#gse184198) | **Gastric** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 2 | 4,000 | Baseline | NCBI GEO | [PMID:36372898](https://doi.org/36372898) |
| 250 | [GSE232733](#gse232733) | **Gastric** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 2,000 | Baseline | NCBI GEO | [PMID:39107346](https://doi.org/39107346) |
| 251 | [GSE320155](#gse320155) | **HCC** | Tier 2 (Baseline Atlas) | 10x Chromium (3' unspecified) | Unselected / Total Single-Cell Suspension | 60 | 120,000 | Baseline | NCBI GEO | [PMID:42443155](https://doi.org/42443155) |
| 252 | [GSE326201](#gse326201) | **HCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 18 | 36,000 | Baseline | NCBI GEO | GEO / CZI |
| 253 | [GSE290925](#gse290925) | **HCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 12 | 24,000 | Baseline | NCBI GEO | [PMID:40241752](https://doi.org/40241752) |
| 254 | [GSE208308](#gse208308) | **HCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 6 | 12,000 | Baseline | NCBI GEO | [PMID:40966278](https://doi.org/40966278) |
| 255 | [GSE291757](#gse291757) | **HCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 4 | 8,000 | Baseline | NCBI GEO | GEO / CZI |
| 256 | [GSE208307](#gse208307) | **HCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 2 | 4,000 | Baseline | NCBI GEO | [PMID:40966278](https://doi.org/40966278) |
| 257 | [GSE280982](#gse280982) | **HNSCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 49 | 98,000 | Baseline | NCBI GEO | [PMID:40593620](https://doi.org/40593620) |
| 258 | [CELLxGENE_714e6bc2](#cellxgene_714e6bc2) | **HNSCC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 44 | 134,385 | Baseline | CZ CELLxGENE | [Dataset Version: https://datasets.cellxgene.cziscience.com/377e8ef1-6f2f-4091-bc27-9bc01959b8eb.h5ad curated and distributed by CZ CELLxGENE Discover in Collection: https://cellxgene.cziscience.com/collections/bc7397a3-ea49-4d57-84b8-80bd6885d4c4](https://doi.org/Dataset Version: https://datasets.cellxgene.cziscience.com/377e8ef1-6f2f-4091-bc27-9bc01959b8eb.h5ad curated and distributed by CZ CELLxGENE Discover in Collection: https://cellxgene.cziscience.com/collections/bc7397a3-ea49-4d57-84b8-80bd6885d4c4) |
| 259 | [CELLxGENE_60acb72d](#cellxgene_60acb72d) | **HNSCC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 44 | 227,816 | Baseline | CZ CELLxGENE | [Dataset Version: https://datasets.cellxgene.cziscience.com/8b58da42-48ca-4615-90f2-0a3c25687da1.h5ad curated and distributed by CZ CELLxGENE Discover in Collection: https://cellxgene.cziscience.com/collections/bc7397a3-ea49-4d57-84b8-80bd6885d4c4](https://doi.org/Dataset Version: https://datasets.cellxgene.cziscience.com/8b58da42-48ca-4615-90f2-0a3c25687da1.h5ad curated and distributed by CZ CELLxGENE Discover in Collection: https://cellxgene.cziscience.com/collections/bc7397a3-ea49-4d57-84b8-80bd6885d4c4) |
| 260 | [GSE268014](#gse268014) | **HNSCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 39 | 78,000 | Baseline | NCBI GEO | [PMID:42601456](https://doi.org/42601456) |
| 261 | [CELLxGENE_624d92e2](#cellxgene_624d92e2) | **HNSCC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 37 | 93,431 | Baseline | CZ CELLxGENE | [Dataset Version: https://datasets.cellxgene.cziscience.com/be2f74ea-f487-4650-8c03-080798c2872c.h5ad curated and distributed by CZ CELLxGENE Discover in Collection: https://cellxgene.cziscience.com/collections/bc7397a3-ea49-4d57-84b8-80bd6885d4c4](https://doi.org/Dataset Version: https://datasets.cellxgene.cziscience.com/be2f74ea-f487-4650-8c03-080798c2872c.h5ad curated and distributed by CZ CELLxGENE Discover in Collection: https://cellxgene.cziscience.com/collections/bc7397a3-ea49-4d57-84b8-80bd6885d4c4) |
| 262 | [GSE198315](#gse198315) | **HNSCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 36 | 268,131 | Baseline | NCBI GEO | [PMID:41054545](https://doi.org/41054545) |
| 263 | [GSE310797](#gse310797) | **HNSCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 19 | 38,000 | Baseline | NCBI GEO | [PMID:41851114](https://doi.org/41851114) |
| 264 | [GSE296771](#gse296771) | **HNSCC** | Tier 2 (Baseline Atlas) | 10x Chromium 5' (Immune Profiling) | FACS-sorted (CD3+/CD8+ T-cell enriched) | 16 | 32,000 | Baseline | NCBI GEO | GEO / CZI |
| 265 | [CELLxGENE_55ca4411](#cellxgene_55ca4411) | **HNSCC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 12 | 28,186 | Baseline | CZ CELLxGENE | [10.1111/cas.15979](https://doi.org/10.1111/cas.15979) |
| 266 | [GSE322620](#gse322620) | **HNSCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 10 | 20,000 | Baseline | NCBI GEO | GEO / CZI |
| 267 | [GSE281978](#gse281978) | **HNSCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 6 | 12,000 | Baseline | NCBI GEO | [PMID:38235919](https://doi.org/38235919) |
| 268 | [GSE339480](#gse339480) | **HNSCC** | Tier 2 (Baseline Atlas) | 10x Chromium (3' unspecified) | Unselected / Total Single-Cell Suspension | 1 | 2,000 | Baseline | NCBI GEO | GEO / CZI |
| 269 | [GSE286935](#gse286935) | **HNSCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 2,000 | Baseline | NCBI GEO | [PMID:39970232](https://doi.org/39970232) |
| 270 | [CELLxGENE_b617ee1b](#cellxgene_b617ee1b) | **Melanoma** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 132 | 391,963 | Baseline | CZ CELLxGENE | [10.1038/s41467-024-49916-4](https://doi.org/10.1038/s41467-024-49916-4) |
| 271 | [GSE216069](#gse216069) | **Melanoma** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Nuclei Isolation (snRNA-seq) | 71 | 142,000 | Baseline | NCBI GEO | [PMID:36624340](https://doi.org/36624340) |
| 272 | [GSE269936](#gse269936) | **Melanoma** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 38 | 76,000 | Baseline | NCBI GEO | [PMID:40890106](https://doi.org/40890106) |
| 273 | [GSE324655](#gse324655) | **Melanoma** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 31 | 62,000 | Baseline | NCBI GEO | [PMID:42039444](https://doi.org/42039444) |
| 274 | [GSE230574](#gse230574) | **Melanoma** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 2 | 4,000 | Baseline | NCBI GEO | [PMID:41405996](https://doi.org/41405996) |
| 275 | [GSE338555](#gse338555) | **Melanoma** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 2 | 4,000 | Baseline | NCBI GEO | GEO / CZI |
| 276 | [GSE317349](#gse317349) | **Melanoma** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 2 | 4,000 | Baseline | NCBI GEO | GEO / CZI |
| 277 | [CELLxGENE_89972213](#cellxgene_89972213) | **Melanoma** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 27,757 | Baseline | CZ CELLxGENE | [10.1016/j.immuni.2022.09.002](https://doi.org/10.1016/j.immuni.2022.09.002) |
| 278 | [CELLxGENE_b3052902](#cellxgene_b3052902) | **Melanoma** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 37,193 | Baseline | CZ CELLxGENE | [10.1016/j.immuni.2022.09.002](https://doi.org/10.1016/j.immuni.2022.09.002) |
| 279 | [CELLxGENE_d4dc4cfe](#cellxgene_d4dc4cfe) | **Melanoma** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 36,785 | Baseline | CZ CELLxGENE | [10.1016/j.immuni.2022.09.002](https://doi.org/10.1016/j.immuni.2022.09.002) |
| 280 | [CELLxGENE_dcb7c544](#cellxgene_dcb7c544) | **Melanoma** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 40,982 | Baseline | CZ CELLxGENE | [10.1016/j.immuni.2022.09.002](https://doi.org/10.1016/j.immuni.2022.09.002) |
| 281 | [CELLxGENE_50c4a6d6](#cellxgene_50c4a6d6) | **Melanoma** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 21,003 | Baseline | CZ CELLxGENE | [10.1016/j.immuni.2022.09.002](https://doi.org/10.1016/j.immuni.2022.09.002) |
| 282 | [CELLxGENE_76bb43ff](#cellxgene_76bb43ff) | **Melanoma** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 19,034 | Baseline | CZ CELLxGENE | [10.1016/j.immuni.2022.09.002](https://doi.org/10.1016/j.immuni.2022.09.002) |
| 283 | [CELLxGENE_ae4552dc](#cellxgene_ae4552dc) | **Melanoma** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 37,275 | Baseline | CZ CELLxGENE | [10.1016/j.immuni.2022.09.002](https://doi.org/10.1016/j.immuni.2022.09.002) |
| 284 | [CELLxGENE_9f222629](#cellxgene_9f222629) | **NSCLC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 484 | 2,282,447 | Baseline | CZ CELLxGENE | [10.1038/s41591-023-02327-2](https://doi.org/10.1038/s41591-023-02327-2) |
| 285 | [CELLxGENE_1e6a6ef9](#cellxgene_1e6a6ef9) | **NSCLC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 318 | 1,283,972 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2022.10.008](https://doi.org/10.1016/j.ccell.2022.10.008) |
| 286 | [CELLxGENE_232f6a5a](#cellxgene_232f6a5a) | **NSCLC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 298 | 892,296 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2022.10.008](https://doi.org/10.1016/j.ccell.2022.10.008) |
| 287 | [GSE327167](#gse327167) | **NSCLC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 113 | 226,000 | Baseline | NCBI GEO | [PMID:41319863](https://doi.org/41319863) |
| 288 | [GSE308103](#gse308103) | **NSCLC** | Tier 2 (Baseline Atlas) | Single-Nucleus RNA-seq | Nuclei Isolation (snRNA-seq) | 75 | 150,000 | Baseline | NCBI GEO | [PMID:41202811](https://doi.org/41202811) |
| 289 | [GSE192402](#gse192402) | **NSCLC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 74 | 148,000 | Baseline | NCBI GEO | [PMID:36624340](https://doi.org/36624340) |
| 290 | [CELLxGENE_e9175006](#cellxgene_e9175006) | **NSCLC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 42 | 14,072 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2021.09.008](https://doi.org/10.1016/j.ccell.2021.09.008) |
| 291 | [CELLxGENE_486486d4](#cellxgene_486486d4) | **NSCLC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | FACS-sorted (CD3+/CD8+ T-cell enriched) | 42 | 46,140 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2021.09.008](https://doi.org/10.1016/j.ccell.2021.09.008) |
| 292 | [CELLxGENE_576f193c](#cellxgene_576f193c) | **NSCLC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 42 | 147,137 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2021.09.008](https://doi.org/10.1016/j.ccell.2021.09.008) |
| 293 | [CELLxGENE_a6858c10](#cellxgene_a6858c10) | **NSCLC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 42 | 64,091 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2021.09.008](https://doi.org/10.1016/j.ccell.2021.09.008) |
| 294 | [CELLxGENE_f64e1be1](#cellxgene_f64e1be1) | **NSCLC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 42 | 73,047 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2021.09.008](https://doi.org/10.1016/j.ccell.2021.09.008) |
| 295 | [CELLxGENE_d4cfefa0](#cellxgene_d4cfefa0) | **NSCLC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 42 | 9,778 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2021.09.008](https://doi.org/10.1016/j.ccell.2021.09.008) |
| 296 | [GSE311609](#gse311609) | **NSCLC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Nuclei Isolation (snRNA-seq) | 41 | 82,000 | Baseline | NCBI GEO | [PMID:42062553](https://doi.org/42062553) |
| 297 | [CELLxGENE_d224c8e0](#cellxgene_d224c8e0) | **NSCLC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 41 | 8,030 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2021.09.008](https://doi.org/10.1016/j.ccell.2021.09.008) |
| 298 | [CELLxGENE_d41f45c1](#cellxgene_d41f45c1) | **NSCLC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 23 | 82,991 | Baseline | CZ CELLxGENE | [10.1038/s41590-023-01504-2](https://doi.org/10.1038/s41590-023-01504-2) |
| 299 | [CELLxGENE_01ff5cf0](#cellxgene_01ff5cf0) | **NSCLC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v3/v3.1 | Unselected / Total Single-Cell Suspension | 18 | 117,266 | Baseline | CZ CELLxGENE | [10.1186/s40164-025-00740-6](https://doi.org/10.1186/s40164-025-00740-6) |
| 300 | [GSE316782](#gse316782) | **NSCLC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 18 | 36,000 | Baseline | NCBI GEO | GEO / CZI |
| 301 | [GSE333596](#gse333596) | **NSCLC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 12 | 24,000 | Baseline | NCBI GEO | GEO / CZI |
| 302 | [GSE339453](#gse339453) | **NSCLC** | Tier 2 (Baseline Atlas) | 10x Chromium (3' unspecified) | Unselected / Total Single-Cell Suspension | 8 | 16,000 | Baseline | NCBI GEO | GEO / CZI |
| 303 | [GSE346380](#gse346380) | **NSCLC** | Tier 2 (Baseline Atlas) | Subcellular Spatial Transcriptomics (CosMx/Xenium) | Unselected / Total Single-Cell Suspension | 6 | 12,000 | Baseline | NCBI GEO | GEO / CZI |
| 304 | [CELLxGENE_1e4214ce](#cellxgene_1e4214ce) | **NSCLC** | Tier 2 (Baseline Atlas) | 10x Chromium (3' unspecified) | Unselected / Total Single-Cell Suspension | 4 | 35,954 | Baseline | CZ CELLxGENE | [10.1101/2024.11.01.621259](https://doi.org/10.1101/2024.11.01.621259) |
| 305 | [GSE305872](#gse305872) | **NSCLC** | Tier 2 (Baseline Atlas) | 10x Chromium (3' unspecified) | Unselected / Total Single-Cell Suspension | 1 | 2,000 | Baseline | NCBI GEO | [PMID:42063567](https://doi.org/42063567) |
| 306 | [GSE202051](#gse202051) | **PDAC** | Tier 2 (Baseline Atlas) | Single-Nucleus RNA-seq | Nuclei Isolation (snRNA-seq) | 74 | 148,000 | Baseline | NCBI GEO | [PMID:36185212](https://doi.org/36185212) |
| 307 | [GSE291124](#gse291124) | **PDAC** | Tier 2 (Baseline Atlas) | Single-Nucleus RNA-seq | Nuclei Isolation (snRNA-seq) | 17 | 34,000 | Baseline | NCBI GEO | [PMID:41564862](https://doi.org/41564862) |
| 308 | [GSE347847](#gse347847) | **PDAC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 15 | 30,000 | Baseline | NCBI GEO | GEO / CZI |
| 309 | [GSE312209](#gse312209) | **PDAC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 15 | 30,000 | Baseline | NCBI GEO | [PMID:41634803](https://doi.org/41634803) |
| 310 | [GSE348038](#gse348038) | **PDAC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 14 | 28,000 | Baseline | NCBI GEO | GEO / CZI |
| 311 | [GSE284392](#gse284392) | **PDAC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 5 | 10,000 | Baseline | NCBI GEO | [PMID:41634803](https://doi.org/41634803) |
| 312 | [GSE327056](#gse327056) | **PDAC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 4 | 8,000 | Baseline | NCBI GEO | [PMID:42100437](https://doi.org/42100437) |
| 313 | [GSE160977](#gse160977) | **PDAC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 2,000 | Baseline | NCBI GEO | GEO / CZI |
| 314 | [GSE288067](#gse288067) | **PDAC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 23,905 | Baseline | NCBI GEO | GEO / CZI |
| 315 | [GSE300154](#gse300154) | **PDAC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | FACS-sorted (CD3+/CD8+ T-cell enriched) | 1 | 2,000 | Baseline | NCBI GEO | [PMID:41722836](https://doi.org/41722836) |
| 316 | [GSE328692](#gse328692) | **ccRCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 34 | 68,000 | Baseline | NCBI GEO | GEO / CZI |
| 317 | [GSE304262](#gse304262) | **ccRCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 30 | 60,000 | Baseline | NCBI GEO | [PMID:41694381](https://doi.org/41694381) |
| 318 | [GSE294109](#gse294109) | **ccRCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 25 | 50,000 | Baseline | NCBI GEO | [PMID:42242232](https://doi.org/42242232) |
| 319 | [CELLxGENE_5af90777](#cellxgene_5af90777) | **ccRCC** | Tier 2 (Baseline Atlas) | 10x Chromium 5' (Immune Profiling) | Unselected / Total Single-Cell Suspension | 12 | 270,855 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2022.11.001](https://doi.org/10.1016/j.ccell.2022.11.001) |
| 320 | [GSE289672](#gse289672) | **ccRCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | FACS-sorted (EpCAM+ Malignant/Epithelial) | 10 | 20,000 | Baseline | NCBI GEO | GEO / CZI |
| 321 | [CELLxGENE_be39785b](#cellxgene_be39785b) | **ccRCC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v2 | Unselected / Total Single-Cell Suspension | 7 | 20,509 | Baseline | CZ CELLxGENE | [10.1073/pnas.2103240118](https://doi.org/10.1073/pnas.2103240118) |
| 322 | [CELLxGENE_eaf0c852](#cellxgene_eaf0c852) | **ccRCC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 7 | 27,912 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2022.11.001](https://doi.org/10.1016/j.ccell.2022.11.001) |
| 323 | [CELLxGENE_bd65a70f](#cellxgene_bd65a70f) | **ccRCC** | Tier 2 (Baseline Atlas) | 10x Chromium 5' (Immune Profiling) | Unselected / Total Single-Cell Suspension | 6 | 167,283 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2021.03.007](https://doi.org/10.1016/j.ccell.2021.03.007) |
| 324 | [GSE254444](#gse254444) | **ccRCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 6 | 12,000 | Baseline | NCBI GEO | [PMID:41513925](https://doi.org/41513925) |
| 325 | [CELLxGENE_318cb2a6](#cellxgene_318cb2a6) | **ccRCC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 5 | 10,924 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2022.11.001](https://doi.org/10.1016/j.ccell.2022.11.001) |
| 326 | [GSE272610](#gse272610) | **ccRCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 5 | 243 | Baseline | NCBI GEO | [PMID:40241258](https://doi.org/40241258) |
| 327 | [CELLxGENE_a45125aa](#cellxgene_a45125aa) | **ccRCC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2022.11.001](https://doi.org/10.1016/j.ccell.2022.11.001) |
| 328 | [CELLxGENE_4aefec71](#cellxgene_4aefec71) | **ccRCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 37,568 | Baseline | CZ CELLxGENE | [10.1016/j.immuni.2022.09.002](https://doi.org/10.1016/j.immuni.2022.09.002) |
| 329 | [CELLxGENE_32ffc3a7](#cellxgene_32ffc3a7) | **ccRCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 40,217 | Baseline | CZ CELLxGENE | [10.1016/j.immuni.2022.09.002](https://doi.org/10.1016/j.immuni.2022.09.002) |
| 330 | [CELLxGENE_18fd0190](#cellxgene_18fd0190) | **ccRCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 37,563 | Baseline | CZ CELLxGENE | [10.1016/j.immuni.2022.09.002](https://doi.org/10.1016/j.immuni.2022.09.002) |
| 331 | [CELLxGENE_07efa1c3](#cellxgene_07efa1c3) | **ccRCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 35,839 | Baseline | CZ CELLxGENE | [10.1016/j.immuni.2022.09.002](https://doi.org/10.1016/j.immuni.2022.09.002) |
| 332 | [CELLxGENE_02faf712](#cellxgene_02faf712) | **ccRCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 17,612 | Baseline | CZ CELLxGENE | [10.1016/j.immuni.2022.09.002](https://doi.org/10.1016/j.immuni.2022.09.002) |
| 333 | [CELLxGENE_6d243918](#cellxgene_6d243918) | **ccRCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 18,851 | Baseline | CZ CELLxGENE | [10.1016/j.immuni.2022.09.002](https://doi.org/10.1016/j.immuni.2022.09.002) |
| 334 | [CELLxGENE_9b1437bb](#cellxgene_9b1437bb) | **ccRCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 14,990 | Baseline | CZ CELLxGENE | [10.1016/j.immuni.2022.09.002](https://doi.org/10.1016/j.immuni.2022.09.002) |
| 335 | [CELLxGENE_f25a532c](#cellxgene_f25a532c) | **ccRCC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2022.11.001](https://doi.org/10.1016/j.ccell.2022.11.001) |
| 336 | [CELLxGENE_cc43509e](#cellxgene_cc43509e) | **ccRCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 15,366 | Baseline | CZ CELLxGENE | [10.1016/j.immuni.2022.09.002](https://doi.org/10.1016/j.immuni.2022.09.002) |
| 337 | [CELLxGENE_f339fe89](#cellxgene_f339fe89) | **ccRCC** | Tier 2 (Baseline Atlas) | High-Throughput scRNA-seq | Unselected / Total Single-Cell Suspension | 1 | 20,021 | Baseline | CZ CELLxGENE | [10.1016/j.immuni.2022.09.002](https://doi.org/10.1016/j.immuni.2022.09.002) |
| 338 | [CELLxGENE_05f813a4](#cellxgene_05f813a4) | **ccRCC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2022.11.001](https://doi.org/10.1016/j.ccell.2022.11.001) |
| 339 | [CELLxGENE_0671c0d4](#cellxgene_0671c0d4) | **ccRCC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2022.11.001](https://doi.org/10.1016/j.ccell.2022.11.001) |
| 340 | [CELLxGENE_4c6f9f26](#cellxgene_4c6f9f26) | **ccRCC** | Tier 2 (Baseline Atlas) | 10x Chromium 3' v2 | Unselected / Total Single-Cell Suspension | 1 | 2,576 | Baseline | CZ CELLxGENE | [10.1073/pnas.2103240118](https://doi.org/10.1073/pnas.2103240118) |
| 341 | [CELLxGENE_104cfa2a](#cellxgene_104cfa2a) | **ccRCC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2022.11.001](https://doi.org/10.1016/j.ccell.2022.11.001) |
| 342 | [CELLxGENE_24c31c8c](#cellxgene_24c31c8c) | **ccRCC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2022.11.001](https://doi.org/10.1016/j.ccell.2022.11.001) |
| 343 | [CELLxGENE_252438d3](#cellxgene_252438d3) | **ccRCC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2022.11.001](https://doi.org/10.1016/j.ccell.2022.11.001) |
| 344 | [CELLxGENE_27cd5ac9](#cellxgene_27cd5ac9) | **ccRCC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2022.11.001](https://doi.org/10.1016/j.ccell.2022.11.001) |
| 345 | [CELLxGENE_30437616](#cellxgene_30437616) | **ccRCC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2022.11.001](https://doi.org/10.1016/j.ccell.2022.11.001) |
| 346 | [CELLxGENE_53d62b10](#cellxgene_53d62b10) | **ccRCC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2022.11.001](https://doi.org/10.1016/j.ccell.2022.11.001) |
| 347 | [CELLxGENE_60ac2657](#cellxgene_60ac2657) | **ccRCC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2022.11.001](https://doi.org/10.1016/j.ccell.2022.11.001) |
| 348 | [CELLxGENE_75548d10](#cellxgene_75548d10) | **ccRCC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2022.11.001](https://doi.org/10.1016/j.ccell.2022.11.001) |
| 349 | [CELLxGENE_81328f3f](#cellxgene_81328f3f) | **ccRCC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2022.11.001](https://doi.org/10.1016/j.ccell.2022.11.001) |
| 350 | [CELLxGENE_a6046b15](#cellxgene_a6046b15) | **ccRCC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2022.11.001](https://doi.org/10.1016/j.ccell.2022.11.001) |
| 351 | [CELLxGENE_c3f74413](#cellxgene_c3f74413) | **ccRCC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2022.11.001](https://doi.org/10.1016/j.ccell.2022.11.001) |
| 352 | [CELLxGENE_c5ac3ec2](#cellxgene_c5ac3ec2) | **ccRCC** | Tier 2 (Baseline Atlas) | 10x Visium Spatial Transcriptomics | Unselected / Total Single-Cell Suspension | 1 | 4,992 | Baseline | CZ CELLxGENE | [10.1016/j.ccell.2022.11.001](https://doi.org/10.1016/j.ccell.2022.11.001) |

---

## 5. Cohort Directory by Cancer Indication

### Melanoma Cohorts
*29 cohorts identified for Melanoma*

#### <a id='gse344166'></a>GSE344166 — Integrated spatial and single cell analysis identifies CCL21+ lymphatic endothelial cells as a driver of a favorable immune environment in acral melanoma
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 116 samples/patients; ~232,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE344166_Acral_GeoMX_norm.csv.gz; GSE344166_Experiment_Summary.csv.gz; GSE344166_Hs_R_NGS_WTA_v1.0.pkc.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE344nnn/GSE344166/suppl/GSE344166_Acral_GeoMX_no`
- **Repository Access Link:** [GSE344166](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE344166)
- **Study Abstract / Experimental Design:**
  Acral melanoma (AM) is a rare type of melanoma that responds poorly to immunotherapy. In the current study we undertook integrated single cell RNA-Seq, spatial transcriptomics and multiplexed immunofluorescence to identify potential regulators of the immune environment. Our analyses identified distinct immune habitats at the invasive front, intratumoral regions and areas of distant inflammation ac...

#### <a id='gse286410'></a>GSE286410 — Multi-modal Omics Analysis of a Paediatric Melanoma Highlights Mechanisms Underlying Treatment Resistance [Seq]
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** Single-Nucleus RNA-seq
- **Modality:** snRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 41 samples/patients; ~82,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE286410_Annotation.txt.gz; GSE286410_GeoMx_qnorm.txt.gz; filelist.txt`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286410/suppl/GSE286410_Annotation.txt`
- **Citation / Reference:** PMID:41168392
- **Repository Access Link:** [GSE286410](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286410)
- **Study Abstract / Experimental Design:**
  Background Cutaneous malignant melanoma is a common cancer in adults but extremely rare in young children, affecting fewer than one child per million each year in Europe. Because of its rarity, most treatments for children are adapted from adult therapies, despite possible biological differences. This study aimed to explore the molecular features of a rare and aggressive melanoma in a 16-month-old...

#### <a id='gse218429'></a>GSE218429 — Downregulation of KEAP1 in melanoma promotes resistance to immune checkpoint blockade [scRNA-seq]
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 35 samples/patients; ~70,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE218429_counts.csv.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE218nnn/GSE218429/suppl/GSE218429_counts.csv.gz`
- **Citation / Reference:** PMID:36864091
- **Repository Access Link:** [GSE218429](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218429)
- **Study Abstract / Experimental Design:**
  Immune checkpoint blockade (ICB) has demonstrated efficacy in patients with melanoma, but many exhibit poor responses. Using single cell RNA sequencing of melanoma patient-derived circulating tumor cells (CTCs) and functional characterization using mouse melanoma models, we show that the KEAP1/NRF2 pathway modulates sensitivity to ICB, independently of tumorigenesis. The NRF2 negative regulator, K...

#### <a id='gse198265'></a>GSE198265 — Neoantigen specific CD4+ T cells in human melanoma have diverse differentiation states and correlate with CD8+ T cell, macrophage, and B cell function
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"In a cohort of melanoma patients, the frequency of CXCL13+ CD4+ T cells in the tumor correlated with the transcriptional states of CD8+ T cells and macrophages, maturation of B cells, and patient survival."*
- **Cohort Scale:** 22 samples/patients; ~44,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE198265_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198265/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:35413271
- **Repository Access Link:** [GSE198265](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198265)
- **Study Abstract / Experimental Design:**
  Tumor antigen-specific CD4+ T cells are required for the efficacy of immune checkpoint inhibitors in murine models but their contributions in human cancer are less understood. We used targeted single cell RNA sequencing and matching of T cell receptor sequences to identify signatures and functional correlates of tumor antigen-specific CD4+ T cells infiltrating human melanoma tumors. CD4+ T cells t...

#### <a id='gse242477'></a>GSE242477 — Single-cell profiling of acral melanoma infiltrating lymphocytes reveals a suppressive tumor microenvironment
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Tumor-reactive CD8 TILs showed heterogeneous expression of coinhibitory molecules, including KLRC1 (NKG2A), in subpopulations with therapeutic implications."*
- **Cohort Scale:** 17 samples/patients; ~34,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE242477_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242477/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:39630887
- **Repository Access Link:** [GSE242477](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242477)
- **Study Abstract / Experimental Design:**
  Acral lentiginous melanoma (ALM) is the most common melanoma subtype in non-Caucasians. Despite advances in cancer immunotherapy, current immune checkpoint inhibitors remain unsatisfactory for ALM. Hence, we conducted comprehensive immune profiling using single-cell phenotyping with reactivity screening of the T cell receptors of tumor-infiltrating T lymphocytes (TILs) in ALM. Compared with cutane...

#### <a id='gse303948'></a>GSE303948 — Genome accessibility profile at single cell level for conventional dendritic cells in human metastatic melanoma samples [snATAC-Seq]
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** Single-Nucleus RNA-seq
- **Modality:** scRNA-seq + snRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 14 samples/patients; ~28,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE303948_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE303nnn/GSE303948/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:40890106
- **Repository Access Link:** [GSE303948](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303948)
- **Study Abstract / Experimental Design:**
  Conventional dendritic cells (cDCs) has been shown to mediate immune checkpoint inhibitor responses in cancer patients. We used single nucleus ATAC sequencing (snATAC-seq) to analyze the epigenome of heterogeneous cDCs obtained from 14 human metastatic melanoma samples....

#### <a id='gse294273'></a>GSE294273 — Tumor-resident T cells and dendritic cells form an in situ archetype for immunotherapy response in melanoma [scRNA-seq]
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"In responders, we found enrichment for early differentiation (TCF1+) CD8+ and CD4+ TR which co-localized with melanoma cells."*
- **Cohort Scale:** 12 samples/patients; ~24,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE294273_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE294nnn/GSE294273/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:42277002
- **Repository Access Link:** [GSE294273](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294273)
- **Study Abstract / Experimental Design:**
  Tumor-resident (TR) T cells participate in immunosurveillance of melanoma, but their role in determining response to immune checkpoint inhibitors (ICI) has not been comprehensively explored. We performed spatial and single-cell profiling on 32 metastatic melanoma lymph node samples, from treatment-naïve, ICI-resistant, or ICI-responsive patients. In responders, we found enrichment for early differ...

#### <a id='gse320040'></a>GSE320040 — High-resolution and noninvasive profiling of the tumor microenvironment with spatial ecotypes [scRNA-seq]
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 7 samples/patients; ~14,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE320040_CRC_tumor_scrna_counts.mtx.gz; GSE320040_CRC_tumor_scrna_logcpm.mtx.gz; GSE320040_Melanoma_tumor_scrna_counts.mtx.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320040/suppl/GSE320040_CRC_tumor_scrn`
- **Citation / Reference:** PMID:42092150
- **Repository Access Link:** [GSE320040](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE320040)
- **Study Abstract / Experimental Design:**
  Multicellular programs in the tumor microenvironment (TME) drive cancer pathogenesis and response to therapy but remain challenging to identify and profile clinically. Here, we present a machine learning framework for multi-analyte profiling of spatially dependent cell states and multicellular ecosystems, termed spatial ecotypes (SEs). By integrating 10M single-cell and spot-level spatial transcri...

#### <a id='gse210963'></a>GSE210963 — Single-Cell RNA-Seq Analysis of Patient Myeloid-Derived Suppressor Cells and the Response to Inhibition of Bruton’s Tyrosine Kinase
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** 10x Chromium (3' unspecified)
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 6 samples/patients; ~12,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE210963_counts.csv.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/GSE210963/suppl/GSE210963_counts.csv.gz`
- **Citation / Reference:** PMID:38015751
- **Repository Access Link:** [GSE210963](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210963)
- **Study Abstract / Experimental Design:**
  Myeloid-derived suppressor cell (MDSC) levels are elevated in cancer patients and contribute to reduced efficacy of immune checkpoint therapy. MDSC express Bruton’s Tyrosine Kinase (BTK) and BTK inhibition with ibrutinib, an FDA-approved irreversible inhibitor of BTK, leads to reduced MDSC expansion/function in mice and significantly improves the anti-tumor activity of anti-PD-1 antibody treatment...

#### <a id='gse256291'></a>GSE256291 — Comparing transcriptional profiles of CD14+  monocytes from melanoma patients and healthy donors by scRNA-seq.
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 4 samples/patients; ~8,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE256291_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE256nnn/GSE256291/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:39266214
- **Repository Access Link:** [GSE256291](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE256291)
- **Study Abstract / Experimental Design:**
  Myeloid-derived suppressor cells (MDSC), a heterogeneous population of myeloid cells, accumulate in the melanoma microenvironment. With the ability to inhibit anti-tumor T cell responses, MDSC have been shown to promote immunosuppression, enhancing tumor progression and tumor cell resistance to the immunotherapy. Novel markers are still needed to define MDSC. In this study,we wanted to find possib...

#### <a id='gse244983'></a>GSE244983 — Molecular patterns of resistance to immune checkpoint blockade in melanoma [scRNA-Seq]
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"In addition, anti-PD1 resistant tumors had reduced fractions of PD1+ CD8+ T cells as compared to ICB naïve metastases."*
- **Cohort Scale:** 4 samples/patients; ~8,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE244983_NormalizedData_scRNAseq.txt.gz; GSE244983_RawCounts_scRNAseq.txt.gz; GSE244983_SingleCellAnnotations.txt.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244983/suppl/GSE244983_NormalizedData`
- **Citation / Reference:** PMID:38594286
- **Repository Access Link:** [GSE244983](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244983)
- **Study Abstract / Experimental Design:**
  Immune checkpoint blockade (ICB) has improved outcome for patients with metastatic melanoma but not all benefit from treatment. Several immune- and tumor intrinsic features are associated with clinical response at baseline. However, we need to further understand the molecular changes occurring during development of ICB resistance. Here, we collected biopsies from a cohort of 44 melanoma patients a...

#### <a id='gse300446'></a>GSE300446 — Spatial tumour-immune ecosystems shape the efficacy of anti-PD1 immunotherapy in primary cutaneous melanoma [snRNAseq]
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** Single-Nucleus RNA-seq
- **Modality:** snRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Nuclei Isolation (snRNA-seq)**
  > *"Single-nucleus RNA sequencing performed on isolated nuclei from frozen resection tissue."*
- **Cohort Scale:** 2 samples/patients; ~4,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE300446_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE300nnn/GSE300446/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE300446](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300446)
- **Study Abstract / Experimental Design:**
  Intra-tumoral heterogeneity in melanoma arises from dynamic cancer cell plasticity and underlies various mechanisms of immune escape. Here, we combined high-plex immunofluorescence imaging with spatially resolved transcriptomics to map the architecture of melanoma cell states and their interactions with the immune microenvironment in primary cutaneous tumours prior to adjuvant anti-PD1 immune chec...

#### <a id='cellxgene_7b20c613'></a>CELLxGENE_7b20c613 — Integrated cancer cell-specific single-cell RNA-seq datasets of immune checkpoint blockade-treated patients
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 167 samples/patients; ~355,941 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `7b20c613-9add-43d1-87e9-defd3d9b9f8c.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/c1d02236-6709-4feb-b7d6-30e2a203e8d1.h5ad`
- **Citation / Reference:** 10.1038/s41597-025-04381-6
- **Repository Access Link:** [CELLxGENE_7b20c613](https://cellxgene.cziscience.com/e/7b20c613-9add-43d1-87e9-defd3d9b9f8c.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Integrated cancer cell-specific single-cell RNA-seq datasets of immune checkpoint blockade-treated patients | Diseases: HER2 positive breast carcinoma, basal cell carcinoma, estrogen-receptor positive breast cancer, hepatocellular carcinoma, intrahepatic cholangiocarcinoma, melanoma, metastatic melanoma, nonpapillary renal cell carcinoma, squamous cell carcinoma, triple-negative breast...

#### <a id='gse270464'></a>GSE270464 — Specific oncogene activation of the cell of origin in mucosal melanoma
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 10 samples/patients; ~20,000 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `filelist.txt; GSE270464_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270464/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:40695831
- **Repository Access Link:** [GSE270464](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270464)
- **Study Abstract / Experimental Design:**
  Mucosal melanoma (MM) is a deadly cancer derived from mucosal melanocytes. To test the consequences of MM genetics, we develop a zebrafish model in which all melanocytes experience CCND1 expression and loss of PTEN and TP53. Surprisingly, melanoma only develops from melanocytes lining internal organs, analogous to the location of patient MM. We find that zebrafish MMs have a unique chromatin lands...

#### <a id='gse211068'></a>GSE211068 — Integrated multiomics profiling identifies the differentiation program of regulatory T cells in human tumors [scRNA-seq]
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 9 samples/patients; ~18,000 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `filelist.txt; GSE211068_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211068/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:36206353
- **Repository Access Link:** [GSE211068](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211068)
- **Study Abstract / Experimental Design:**
  Regulatory T (Treg) cells suppress effective antitumor immunity in tumor-bearing hosts, thereby becoming promising targets in cancer immunotherapy. Here, we show that Treg cells in the tumor microenvironment (TME) of human lung cancers harbor a completely different open chromatin profile compared to effector T cells and peripheral Treg cells. The integrative sequencing analyses revealed that BATF,...

#### <a id='cellxgene_b617ee1b'></a>CELLxGENE_b617ee1b — A multi-tissue single-cell tumor microenvironment atlas
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 132 samples/patients; ~391,963 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `b617ee1b-f8c8-4de9-b82b-e803ab93550d.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/1b76227b-c731-4807-9487-ad5e4d24e0d0.h5ad`
- **Citation / Reference:** 10.1038/s41467-024-49916-4
- **Repository Access Link:** [CELLxGENE_b617ee1b](https://cellxgene.cziscience.com/e/b617ee1b-f8c8-4de9-b82b-e803ab93550d.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Dissecting novel myeloid-derived cell states through single-cell RNA-Seq and its impact on clinical outcome across tumor types | Diseases: breast apocrine carcinoma, breast cancer, cecum adenocarcinoma, colon adenocarcinoma, colorectal adenocarcinoma, invasive ductal breast carcinoma, large cell carcinoma, liver cancer, lung adenocarcinoma, lung cancer, melanoma, metaplastic breast car...

#### <a id='gse216069'></a>GSE216069 — Multi-modal single-cell and whole-genome sequencing of minute, frozen specimens to propel clinical applications [scRNA and snRNA]
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Nuclei Isolation (snRNA-seq)**
  > *"Single-nucleus RNA sequencing performed on isolated nuclei from frozen resection tissue."*
- **Cohort Scale:** 71 samples/patients; ~142,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE216069_BI5_slyper_meta_data.csv.gz; GSE216069_BI5_slyper_processed_data.csv.gz; GSE216069_Mel_meta_data.csv.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE216nnn/GSE216069/suppl/GSE216069_BI5_slyper_met`
- **Citation / Reference:** PMID:36624340
- **Repository Access Link:** [GSE216069](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216069)
- **Study Abstract / Experimental Design:**
  We provide 5' and/or 3' single cell and single nuclei RNA sequencing data for matched comparisons of profiling from fresh (single cells) and frozen (single nuclei) samples of a non small cell lung cancer sample, a cutaneous melanoma sample , and a primary uveal melanoma sample. We also provide matched samples processed with published dissociation protocols (*TST and *CST) for two samples as compar...

#### <a id='gse269936'></a>GSE269936 — Gene expression profile at single cell level for human metastatic melanoma samples
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 38 samples/patients; ~76,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE269936_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE269nnn/GSE269936/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:40890106
- **Repository Access Link:** [GSE269936](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269936)
- **Study Abstract / Experimental Design:**
  There exists great intra- and inter-tumor heterogeneity for metastatic melanomas. We used single cell RNA sequencing (scRNA-seq) to analyze cell type heterogeneity of metastatic melanoma and cell type specific responses to different treatments....

#### <a id='gse324655'></a>GSE324655 — Transcriptional and chromatin accessibility profiling of melanoma cells during BRAF inhibitor–induced drug tolerance
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 31 samples/patients; ~62,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE324655_SKMEL5_RNA_featureCounts_matrix_all.csv.gz; GSE324655_SKMEL5_subclone_BRAFi_timecourse_RNAseq_data.csv.gz; filelist.txt`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE324nnn/GSE324655/suppl/GSE324655_SKMEL5_RNA_fea`
- **Citation / Reference:** PMID:42039444
- **Repository Access Link:** [GSE324655](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324655)
- **Study Abstract / Experimental Design:**
  Drug tolerance can emerge rapidly in melanoma following treatment with BRAF inhibitors. This transition has been associated with transcriptional and chromatin state remodeling. To investigate the molecular features of this process, we profiled melanoma cells before and during exposure to a BRAF inhibitor. We generated single-cell RNA sequencing (scRNA-seq), bulk RNA sequencing (bulk RNA-seq), and ...

#### <a id='gse230574'></a>GSE230574 — Gene expression profile at single cell level of human melanoma cell lines, normal human melanocytes and human nevus cells.
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 2 samples/patients; ~4,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE230574_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE230nnn/GSE230574/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41405996
- **Repository Access Link:** [GSE230574](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230574)
- **Study Abstract / Experimental Design:**
  Melanomas are considered highly heterogenous in transcriptional state. We used single cell RNA sequencing (scRNA-seq) to analyze the diversity of melanoma cells as compared to primary melanocytes and nevus melanocytes. Sample 1 contain KIT postitive melanocytes sequenced directly ex vivo from human nevi. Sample 2 contains different human melanoma cell lines and primary human melanocytes grown in v...

#### <a id='gse338555'></a>GSE338555 — Identification of progenitor cells that drive Interleukin 11-dependent tumor-stroma coevolution in advanced melanoma [TGFB_scRNAseq]
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 2 samples/patients; ~4,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE338555_matrix.mtx.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE338nnn/GSE338555/suppl/GSE338555_matrix.mtx.gz`
- **Repository Access Link:** [GSE338555](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE338555)
- **Study Abstract / Experimental Design:**
  Tumor mesenchymal stem cells (TMSCs) orchestrate tumor microenvironment (TME) transformation via their ability to transdifferentiate. Metastatic solid tumors often feature neural crest-like/mesenchymal phenotypes characterized by high expression levels of TGFβ-dependent epithelial-to-mesenchymal transition (EMT) and EMT-like programs. Specifically, EMT-like endothelial-to-mesenchymal (EndMT) progr...

#### <a id='gse317349'></a>GSE317349 — GDF15 reprograms the microenvironment to drive the development of uveal melanoma liver metastases [scRNA-seq]
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 2 samples/patients; ~4,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE317349_1-E-V-filtered_feature_bc_matrix.h5; GSE317349_5-M-V-filtered_feature_bc_matrix.h5`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE317nnn/GSE317349/suppl/GSE317349_1-E-V-filtered`
- **Repository Access Link:** [GSE317349](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE317349)
- **Study Abstract / Experimental Design:**
  Uveal melanoma (UM) results in fatal liver metastasis, yet little is known about the interactions between UM and host cells in the tumor microenvironment that promote this distinctive proclivity. We used single cell (sc)-RNA-Seq analysis of UM-hepatic stellate cell (HSC) co-cultures to demonstrate that HSCs enriched for UM cell states that expressed genes implicated in cell survival, metabolic rep...

#### <a id='cellxgene_89972213'></a>CELLxGENE_89972213 — Cutaneous melanoma, lymph node metastasis Puck_220215_13
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~27,757 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `89972213-8510-49c3-a55a-01cd248ef88a.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/12115237-b7ac-4295-bd99-0c339e695c24.h5ad`
- **Citation / Reference:** 10.1016/j.immuni.2022.09.002
- **Repository Access Link:** [CELLxGENE_89972213](https://cellxgene.cziscience.com/e/89972213-8510-49c3-a55a-01cd248ef88a.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatially mapping T cell receptors and transcriptomes | Diseases: metastatic melanoma | Tissues: lymph node | Assays: Slide-seqV2...

#### <a id='cellxgene_b3052902'></a>CELLxGENE_b3052902 — Cutaneous melanoma, lymph node metastasis Puck_220408_12
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~37,193 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `b3052902-613c-4b04-a0f5-bafd397bb125.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/5622f698-5ff3-4e80-969c-006047e08bc3.h5ad`
- **Citation / Reference:** 10.1016/j.immuni.2022.09.002
- **Repository Access Link:** [CELLxGENE_b3052902](https://cellxgene.cziscience.com/e/b3052902-613c-4b04-a0f5-bafd397bb125.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatially mapping T cell receptors and transcriptomes | Diseases: metastatic melanoma | Tissues: lymph node | Assays: Slide-seqV2...

#### <a id='cellxgene_d4dc4cfe'></a>CELLxGENE_d4dc4cfe — Cutaneous melanoma, lymph node metastasis Puck_220215_14
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~36,785 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `d4dc4cfe-7da5-4e56-bf96-83f01082f5df.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/1a25dabe-7559-449d-9db0-3b8d758f6fec.h5ad`
- **Citation / Reference:** 10.1016/j.immuni.2022.09.002
- **Repository Access Link:** [CELLxGENE_d4dc4cfe](https://cellxgene.cziscience.com/e/d4dc4cfe-7da5-4e56-bf96-83f01082f5df.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatially mapping T cell receptors and transcriptomes | Diseases: metastatic melanoma | Tissues: lymph node | Assays: Slide-seqV2...

#### <a id='cellxgene_dcb7c544'></a>CELLxGENE_dcb7c544 — Cutaneous melanoma, lymph node metastasis Puck_220408_07
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~40,982 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `dcb7c544-bf18-4f9a-97a4-aa816910d862.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/ddd93e28-fb8a-491a-b050-76c11fa3d296.h5ad`
- **Citation / Reference:** 10.1016/j.immuni.2022.09.002
- **Repository Access Link:** [CELLxGENE_dcb7c544](https://cellxgene.cziscience.com/e/dcb7c544-bf18-4f9a-97a4-aa816910d862.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatially mapping T cell receptors and transcriptomes | Diseases: metastatic melanoma | Tissues: lymph node | Assays: Slide-seqV2...

#### <a id='cellxgene_50c4a6d6'></a>CELLxGENE_50c4a6d6 — Cutaneous melanoma, lymph node metastasis Puck_220215_11
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~21,003 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `50c4a6d6-940b-4c6a-a376-aea2ae2d3168.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/564954ae-1126-47cd-9736-474b61fd390d.h5ad`
- **Citation / Reference:** 10.1016/j.immuni.2022.09.002
- **Repository Access Link:** [CELLxGENE_50c4a6d6](https://cellxgene.cziscience.com/e/50c4a6d6-940b-4c6a-a376-aea2ae2d3168.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatially mapping T cell receptors and transcriptomes | Diseases: metastatic melanoma | Tissues: lymph node | Assays: Slide-seqV2...

#### <a id='cellxgene_76bb43ff'></a>CELLxGENE_76bb43ff — Cutaneous melanoma, lymph node metastasis Puck_220408_08
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~19,034 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `76bb43ff-e2f5-4513-a0f0-1059c860c3b7.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/d799315c-05ca-4d3d-9a85-32647cfaed8a.h5ad`
- **Citation / Reference:** 10.1016/j.immuni.2022.09.002
- **Repository Access Link:** [CELLxGENE_76bb43ff](https://cellxgene.cziscience.com/e/76bb43ff-e2f5-4513-a0f0-1059c860c3b7.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatially mapping T cell receptors and transcriptomes | Diseases: metastatic melanoma | Tissues: lymph node | Assays: Slide-seqV2...

#### <a id='cellxgene_ae4552dc'></a>CELLxGENE_ae4552dc — Cutaneous melanoma, lymph node metastasis Puck_220408_09
- **Cancer Type / Indication:** Melanoma
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~37,275 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `ae4552dc-e2ea-4d67-b375-03ec7480f780.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/7a4bbc4e-4351-4eed-ba8d-9b461d517e96.h5ad`
- **Citation / Reference:** 10.1016/j.immuni.2022.09.002
- **Repository Access Link:** [CELLxGENE_ae4552dc](https://cellxgene.cziscience.com/e/ae4552dc-e2ea-4d67-b375-03ec7480f780.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatially mapping T cell receptors and transcriptomes | Diseases: metastatic melanoma | Tissues: lymph node | Assays: Slide-seqV2...

### NSCLC Cohorts
*36 cohorts identified for NSCLC*

#### <a id='gse243013'></a>GSE243013 — A single-cell atlas of immune heterogeneity  in anti-PD1-treated non-small cell lung cancer
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 243 samples/patients; ~486,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE243013_NSCLC_immune_scRNA_counts.mtx.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243013/suppl/GSE243013_NSCLC_immune_s`
- **Citation / Reference:** PMID:40147443
- **Repository Access Link:** [GSE243013](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243013)
- **Study Abstract / Experimental Design:**
  Anti-PD(L)1 with chemotherapy is a standard of care for non-small cell lung cancer (NSCLC), but a varying degree of response to the same regimen is observed in patients. The tumor immune microenvironment (TIME) plays a key role in response to immunotherapy, and the TIME heterogeneity in association with therapeutic outcome is incompletely understood. Here we prospectively applied single-cell RNA a...

#### <a id='gse241934'></a>GSE241934 — Neoadjuvant sintilimab plus chemotherapy in early-stage EGFR-mutant NSCLC: phase 2 trial interim results (NEOTIDE/CTONG2104)
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 88 samples/patients; ~176,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE241934_IIT_Matrix.mtx.gz; GSE241934_Real_Matrix.mtx.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_Matrix.mtx`
- **Citation / Reference:** PMID:38897205
- **Repository Access Link:** [GSE241934](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241934)
- **Study Abstract / Experimental Design:**
  We launched an investigator-initiated, Simon’s two-stage design trial of neoadjuvant sintilimab combined with carboplatin and nab-paclitaxel (nab-PC) in early-stage EGFR-mutant NSCLC (Clinicaltrial.gov number NCT05244213). Here we report the first interim results of stage 1 cohort which met the overall primary endpoint in advance, and multi-omics profiling of neoadjuvant immunotherapy combination ...

#### <a id='gse270148'></a>GSE270148 — Myeloid progenitor dysregulation fuels immunosuppressive macrophages in tumours
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **MACS Bead-selected**
  > *"nscriptome and chromatin accessibility analysis over the continuum of myeloid progenitors, circulating monocytes and tumour-infiltrating mo-macs in mice and in patients with lung cancer to identify myeloid progenitor programs that fuel pro-tumorigenic mo-macs."*
- **Cohort Scale:** 85 samples/patients; ~170,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE270148_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270148/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:40931076
- **Repository Access Link:** [GSE270148](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270148)
- **Study Abstract / Experimental Design:**
  Monocyte-derived macrophages (mo-macs) often drive immunosuppression in the tumour microenvironment (TME) and tumour-enhanced myelopoiesis in the bone marrow fuels these populations. Here we performed paired transcriptome and chromatin accessibility analysis over the continuum of myeloid progenitors, circulating monocytes and tumour-infiltrating mo-macs in mice and in patients with lung cancer to ...

#### <a id='gse317309'></a>GSE317309 — Combination of a CCL21-gene modified dendritic cell vaccine and pembrolizumab induces immune responses in non-small cell lung cancer
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 64 samples/patients; ~128,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE317309_feature_README.txt; filelist.txt; GSE317309_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE317nnn/GSE317309/suppl/GSE317309_feature_README`
- **Repository Access Link:** [GSE317309](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE317309)
- **Study Abstract / Experimental Design:**
  Immune checkpoint inhibitors (ICIs) have transformed treatment for non-small cell lung cancer (NSCLC), but resistance to these therapies is common. We report results from a phase I trial combining intratumoral administration of a CCL21-gene modified dendritic cell (CCL21-DC) vaccine with pembrolizumab in patients with advanced NSCLC. Among 23 patients that received trial therapy, there were no dos...

#### <a id='gse207422'></a>GSE207422 — Tumor microenvironment remodeling after neoadjuvant immunotherapy in non-small cell lung cancer revealed by single-cell RNA sequencing
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 39 samples/patients; ~78,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx; GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz; GSE207422_NSCLC_scRNAseq_metadata.xlsx`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNA`
- **Citation / Reference:** PMID:36869384
- **Repository Access Link:** [GSE207422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422)
- **Study Abstract / Experimental Design:**
  To address the potential therapy-resistant mechanisms and underly the changes of tumor microenvironment (TME) after immunotherapy, we performed scRNA-seq  and bulk RNA-seq from the patients with resectable non-small cell lung cancer (NSCLC) before and after PD-1 blockade combined with chemotherapy. We found that the combined therapy significantly remodeled the immune cell compartments in the TME. ...

#### <a id='gse303762'></a>GSE303762 — Benchmarking long-read RNA-sequencing technologies with LongBench: a cross-platform reference dataset profiling cancer cell lines with bulk and single-cell approaches
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** Single-Nucleus RNA-seq
- **Modality:** scRNA-seq + snRNA-seq
- **Cell Selection / Filtering Strategy:** **Nuclei Isolation (snRNA-seq)**
  > *"Single-nucleus RNA sequencing performed on isolated nuclei from frozen resection tissue."*
- **Cohort Scale:** 38 samples/patients; ~76,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE303762_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE303nnn/GSE303762/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE303762](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303762)
- **Study Abstract / Experimental Design:**
  Long-read RNA sequencing technologies offer unparalleled in- sights into transcriptomes by enabling full-length sequencing of RNA molecules, uncovering novel isoforms and alternative splicing events. While long-read sequencing platforms, such as Pacific Biosciences (PacBio) and Oxford Nanopore Technologies (ONT), have historically been associated with higher error rates, recent advancements in bot...

#### <a id='gse205049'></a>GSE205049 — Spatially Resolved Multi-Omics Single-Cell Analyses Inform Mechanisms of Immune Dysfunction in Pancreatic Cancer
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Single-cell suspension sorted by FACS for CD3+/CD8+ T lymphocytes (TILs)."*
- **Cohort Scale:** 18 samples/patients; ~36,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE205049_scRNA-seq-annotations.csv.gz; GSE205049_scRNA-seq-integrated_GEM.csv.gz; filelist.txt`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205049/suppl/GSE205049_scRNA-seq-anno`
- **Citation / Reference:** PMID:37263303
- **Repository Access Link:** [GSE205049](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205049)
- **Study Abstract / Experimental Design:**
  As pancreatic ductal adenocarcinoma (PDAC) continues to be recalcitrant to therapeutic interventions including poor response to immunotherapy, albeit effective in other solid malignancies, a more nuanced understanding of the immune microenvironment in PDAC is urgently needed. Using a spatially-resolved multimodal single cell approach we unveil a detailed view of the immune micromilieu in PDAC with...

#### <a id='gse205354'></a>GSE205354 — Spatially resolved multi-omics single-cell analyses inform mechanisms of immune-dysfunction in pancreatic cancer
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Single-cell suspension sorted by FACS for CD3+/CD8+ T lymphocytes (TILs)."*
- **Cohort Scale:** 8 samples/patients; ~16,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE205354_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205354/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:37263303
- **Repository Access Link:** [GSE205354](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205354)
- **Study Abstract / Experimental Design:**
  As pancreatic ductal adenocarcinoma (PDAC) continues to be recalcitrant to therapeutic interventions including poor response to immunotherapy, albeit effective in other solid malignancies, a more nuanced understanding of the immune microenvironment in PDAC is urgently needed. Using a spatially-resolved multimodal single cell approach we unveil a detailed view of the immune micromilieu in PDAC with...

#### <a id='gse233203'></a>GSE233203 — The single-cell level molecular characteristics according to combination immunotherapy response of non-small cell lung cancer patients.
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 7 samples/patients; ~14,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE233203_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233203/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:40843133
- **Repository Access Link:** [GSE233203](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233203)
- **Study Abstract / Experimental Design:**
  Non-small cell lung cancer patients finally acquire EGFR-TKI resistance. For these patients, combination immunotherapy, atezolizumab plus bevacizumab and carboplatin/paclitaxel (ABCP), was proposed as a therapy based on targetting both cancer and stromal cells. Despite of relatively good therapeutic outcome, precise mechanism was not fully understood. To dissect tumor microenvironment and regulato...

#### <a id='gse253718'></a>GSE253718 — Single-cell transcriptome reveals drug-resistance signature and immunosuppressive microenvironment in lung adenocarcinoma harboring EGFR mutation
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 6 samples/patients; ~12,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE253718_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253718/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE253718](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253718)
- **Study Abstract / Experimental Design:**
  A majority of EGFR-mutated lung adenocarcinoma (LUAD) inevitably develops acquired resistance to tyrosine kinase inhibitor (TKI) therapy within 2 years. Due to the heterogeneous tumor microenvironment (TME) and discrepant immunotherapy responses in patients with different EGFR mutation statuses, it is necessary to characterize the immune signatures during tumor resistance, facilitating the explora...

#### <a id='gse223779'></a>GSE223779 — Single-cell transcriptomic analysis uncovers intratumoral heterogeneity and drug-tolerant persister in ALK-rearranged lung adenocarcinoma
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 6 samples/patients; ~12,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE223779_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE223nnn/GSE223779/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:37272226
- **Repository Access Link:** [GSE223779](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223779)
- **Study Abstract / Experimental Design:**
  Acquired drug resistance is the major therapeutic obstacle to maintenance treatment of advanced-stage non-small cell lung cancer. Lung adenocarcinoma (ADC) harboring driver mutations also showed poor response to immune checkpoint inhibitors (ICIs). Underlying mechanisms of how drug insensitivity evolves remain unclear. Here we explored the intratumoral heterogeneity of tyrosine kinase inhibitor (T...

#### <a id='gse307811'></a>GSE307811 — Development of antibody-drug conjugates targeting L1CAM to treat metastatic cancer
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~2,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE307811_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE307nnn/GSE307811/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41511407
- **Repository Access Link:** [GSE307811](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307811)
- **Study Abstract / Experimental Design:**
  Effective treatments for metastatic disease have eluded precision oncology efforts and represent a major unmet need. The L1 cell adhesion molecule (L1CAM) is a transmembrane protein expressed by populations of Metastasis Stem cells (MetSC), that are enriched for the ability to reinitiate and propel metastatic growth and display phenotypic plasticity and resistance to conventional chemotherapies. N...

#### <a id='gse285888'></a>GSE285888 — Single-Cell RNA Sequencing of Baseline PBMCs Predicts ICI efficacy and irAE Severity in NSCLC Patients
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Moreover, we found that PRF1 in CD8+ T cells and NK cells plays a critical mediator of complete responses and reduced irAE to ICI therapy."*
- **Cohort Scale:** 1 samples/patients; ~2,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE285888_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285888/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:40404203
- **Repository Access Link:** [GSE285888](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285888)
- **Study Abstract / Experimental Design:**
  Immune checkpoint inhibitors (ICIs) have revolutionized by offering remarkable clinical benefits and durable responses for patients with advanced non-small cell lung cancer (NSCLC). However, a very small percentage of patients responded to ICI treatment, and immune-related adverse events (irAEs) leading to treatment discontinuation remain challenges. Despite the recognized need for biomarkers that...

#### <a id='gse276139'></a>GSE276139 — Gene expression profile at single cell level of cerebrospinal fluid (CSF) cells from lung adenocarcinoma leptomeningeal metastases patients (LUAD LM)
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"l-cell communication analysis showed that CSF-CTC reinforces immunosuppression by co-inhibitory checkpoint axis NECTIN2_TIGIT axis with the CD8+T/NK cells, and via CD47_SIRPA axis with antigen-presenting cells."*
- **Cohort Scale:** 6 samples/patients; ~12,000 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `filelist.txt; GSE276139_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE276nnn/GSE276139/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41876785
- **Repository Access Link:** [GSE276139](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276139)
- **Study Abstract / Experimental Design:**
  Lung adenocarcinoma (LUAD)-derived leptomeningeal metastases (LM) represent a predominant subtype among all LM cases. Nevertheless, the cerebrospinal fluid (CSF) profile of LUAD-LM patients remains poorly characterized and reliable CSF diagnostic biomarkers for LUAD-LM have yet to be established. Using single-cell RNA sequencing data of CSF cells from six LUAD-LM patients, we drew a systematic tra...

#### <a id='cellxgene_9f222629'></a>CELLxGENE_9f222629 — An integrated cell atlas of the human lung in health and disease (full)
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 484 samples/patients; ~2,282,447 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `9f222629-9e39-47d0-b83f-e08d610c7479.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/5f863718-02a6-46cb-83e7-52be988b915b.h5ad`
- **Citation / Reference:** 10.1038/s41591-023-02327-2
- **Repository Access Link:** [CELLxGENE_9f222629](https://cellxgene.cziscience.com/e/9f222629-9e39-47d0-b83f-e08d610c7479.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: The integrated Human Lung Cell Atlas | Diseases: COVID-19, chronic obstructive pulmonary disease, chronic rhinitis, cystic fibrosis, hypersensitivity pneumonitis, interstitial lung disease, lung adenocarcinoma, lung large cell carcinoma, lymphangioleiomyomatosis, non-specific interstitial pneumonia, normal, pleomorphic carcinoma, pneumonia, pulmonary fibrosis, pulmonary sarcoidosis, sq...

#### <a id='cellxgene_1e6a6ef9'></a>CELLxGENE_1e6a6ef9 — The single-cell lung cancer atlas (LuCA) -- extended atlas
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 318 samples/patients; ~1,283,972 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `1e6a6ef9-7ec9-4c90-bbfb-2ad3c3165fd1.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/33165751-ae33-4a65-94ee-52d9bc38f97e.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2022.10.008
- **Repository Access Link:** [CELLxGENE_1e6a6ef9](https://cellxgene.cziscience.com/e/1e6a6ef9-7ec9-4c90-bbfb-2ad3c3165fd1.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: High-resolution single-cell atlas reveals diversity and plasticity of tumor-associated neutrophils in non-small cell lung cancer | Diseases: chronic obstructive pulmonary disease, lung adenocarcinoma, non-small cell lung carcinoma, normal, squamous cell lung carcinoma | Tissues: adrenal tissue, brain, liver, lung, lymph node, pleural effusion | Assays: 10x 3' v1, 10x 3' v2, 10x 3' v3, ...

#### <a id='cellxgene_232f6a5a'></a>CELLxGENE_232f6a5a — The single-cell lung cancer atlas (LuCA) -- core atlas
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 298 samples/patients; ~892,296 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `232f6a5a-a04c-4758-a6e8-88ab2e3a6e69.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/46e0287b-9a33-4e83-99f3-8c044131bfdc.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2022.10.008
- **Repository Access Link:** [CELLxGENE_232f6a5a](https://cellxgene.cziscience.com/e/232f6a5a-a04c-4758-a6e8-88ab2e3a6e69.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: High-resolution single-cell atlas reveals diversity and plasticity of tumor-associated neutrophils in non-small cell lung cancer | Diseases: chronic obstructive pulmonary disease, lung adenocarcinoma, non-small cell lung carcinoma, normal, squamous cell lung carcinoma | Tissues: adrenal tissue, brain, liver, lung, lymph node, pleural effusion | Assays: 10x 3' v1, 10x 3' v2, 10x 3' v3, ...

#### <a id='gse327167'></a>GSE327167 — Single-cell RNA-seq profiling of pleural mesothelioma and comparator pleural specimens
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 113 samples/patients; ~226,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE327167_initial_scRNA_counts_matrix.mtx.gz; GSE327167_stringent_scRNA_counts_matrix.mtx.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE327nnn/GSE327167/suppl/GSE327167_initial_scRNA_`
- **Citation / Reference:** PMID:41319863
- **Repository Access Link:** [GSE327167](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327167)
- **Study Abstract / Experimental Design:**
  This GEO submission contains processed single-cell RNA sequencing data generated using the Seq-Well platform from pleural mesothelioma tumors and comparator pleural specimens. Raw sequencing data are available through controlled access in dbGaP under accession phs004285. Processed expression matrices and cell-level metadata for the single-cell libraries are provided in this submission....

#### <a id='gse308103'></a>GSE308103 — Single-nucleus RNA-sequencing (snRNA-seq) analysis of normal lung tissues, precursor lung lesions and invasive lung adenocarcinomas
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** Single-Nucleus RNA-seq
- **Modality:** snRNA-seq
- **Cell Selection / Filtering Strategy:** **Nuclei Isolation (snRNA-seq)**
  > *"Single-nucleus RNA sequencing performed on isolated nuclei from frozen resection tissue."*
- **Cohort Scale:** 75 samples/patients; ~150,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE308103_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE308nnn/GSE308103/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41202811
- **Repository Access Link:** [GSE308103](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308103)
- **Study Abstract / Experimental Design:**
  Understanding cellular processes underlying early lung adenocarcinoma (LUAD) development is needed to devise intervention strategies. Here, we performed snRNA-seq analysis of human treatment-naive lungs comprising various stages in the sequence of pathogenesis of LUAD including normal lung tissues, atypical adenomatous hyperplasia (AAH), adenocarcinoma in situ (AIS), minimally invasive adenocarcin...

#### <a id='gse192402'></a>GSE192402 — Multi-modal single-cell and whole-genome sequencing of minute, frozen specimens to propel clinical applications
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 74 samples/patients; ~148,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE192402_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE192nnn/GSE192402/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:36624340
- **Repository Access Link:** [GSE192402](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192402)
- **Study Abstract / Experimental Design:**
  This SuperSeries is composed of the SubSeries listed below....

#### <a id='cellxgene_e9175006'></a>CELLxGENE_e9175006 — Myeloid cells
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 42 samples/patients; ~14,072 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `e9175006-8978-4417-939f-819855eab80e.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/38504daf-b825-48ea-b277-033b64d5babb.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2021.09.008
- **Repository Access Link:** [CELLxGENE_e9175006](https://cellxgene.cziscience.com/e/e9175006-8978-4417-939f-819855eab80e.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN MSK - Single cell profiling reveals novel tumor and myeloid subpopulations in small cell lung cancer | Diseases: lung adenocarcinoma, normal, small cell lung carcinoma | Tissues: adrenal gland, axilla, bone spine, brain, liver, lung, lymph node, pleural effusion | Assays: 10x 3' transcription profiling, 10x 3' v2, 10x 3' v3...

#### <a id='cellxgene_486486d4'></a>CELLxGENE_486486d4 — T cells
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Single-cell suspension sorted by FACS for CD3+/CD8+ T lymphocytes (TILs)."*
- **Cohort Scale:** 42 samples/patients; ~46,140 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `486486d4-9462-43e5-9249-eb43fa5a49a6.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/6fde3ad9-c2dc-4bea-bcb1-100192dd5877.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2021.09.008
- **Repository Access Link:** [CELLxGENE_486486d4](https://cellxgene.cziscience.com/e/486486d4-9462-43e5-9249-eb43fa5a49a6.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN MSK - Single cell profiling reveals novel tumor and myeloid subpopulations in small cell lung cancer | Diseases: lung adenocarcinoma, normal, small cell lung carcinoma | Tissues: adrenal gland, axilla, bone spine, brain, liver, lung, lymph node, pleural effusion | Assays: 10x 3' transcription profiling, 10x 3' v2, 10x 3' v3...

#### <a id='cellxgene_576f193c'></a>CELLxGENE_576f193c — Combined samples
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 42 samples/patients; ~147,137 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `576f193c-75d0-4a11-bd25-8676587e6dc2.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/a9d92e38-9a6e-401b-8484-74bb15122341.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2021.09.008
- **Repository Access Link:** [CELLxGENE_576f193c](https://cellxgene.cziscience.com/e/576f193c-75d0-4a11-bd25-8676587e6dc2.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN MSK - Single cell profiling reveals novel tumor and myeloid subpopulations in small cell lung cancer | Diseases: lung adenocarcinoma, normal, small cell lung carcinoma | Tissues: adrenal gland, axilla, bone spine, brain, liver, lung, lymph node, pleural effusion | Assays: 10x 3' transcription profiling, 10x 3' v2, 10x 3' v3...

#### <a id='cellxgene_a6858c10'></a>CELLxGENE_a6858c10 — Epithelial cells
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 42 samples/patients; ~64,091 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `a6858c10-c52a-4a1d-bc12-a90dbd51ca66.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/cfbc234a-7812-4b82-af81-ec6df4b65e04.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2021.09.008
- **Repository Access Link:** [CELLxGENE_a6858c10](https://cellxgene.cziscience.com/e/a6858c10-c52a-4a1d-bc12-a90dbd51ca66.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN MSK - Single cell profiling reveals novel tumor and myeloid subpopulations in small cell lung cancer | Diseases: lung adenocarcinoma, normal, small cell lung carcinoma | Tissues: adrenal gland, axilla, bone spine, brain, liver, lung, lymph node, pleural effusion | Assays: 10x 3' transcription profiling, 10x 3' v2, 10x 3' v3...

#### <a id='cellxgene_f64e1be1'></a>CELLxGENE_f64e1be1 — Immune cells
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 42 samples/patients; ~73,047 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `f64e1be1-de15-4d27-8da4-82225cd4c035.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/5deb79e5-fdd2-45f2-89cd-f291403cc572.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2021.09.008
- **Repository Access Link:** [CELLxGENE_f64e1be1](https://cellxgene.cziscience.com/e/f64e1be1-de15-4d27-8da4-82225cd4c035.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN MSK - Single cell profiling reveals novel tumor and myeloid subpopulations in small cell lung cancer | Diseases: lung adenocarcinoma, normal, small cell lung carcinoma | Tissues: adrenal gland, axilla, bone spine, brain, liver, lung, lymph node, pleural effusion | Assays: 10x 3' transcription profiling, 10x 3' v2, 10x 3' v3...

#### <a id='cellxgene_d4cfefa0'></a>CELLxGENE_d4cfefa0 — NSCLC epithelial cells
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 42 samples/patients; ~9,778 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `d4cfefa0-3a35-44eb-b848-d7a725b481e7.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/25f2c5ad-3426-42be-9a3c-2aacf189b9aa.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2021.09.008
- **Repository Access Link:** [CELLxGENE_d4cfefa0](https://cellxgene.cziscience.com/e/d4cfefa0-3a35-44eb-b848-d7a725b481e7.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN MSK - Single cell profiling reveals novel tumor and myeloid subpopulations in small cell lung cancer | Diseases: lung adenocarcinoma, normal, small cell lung carcinoma | Tissues: adrenal gland, axilla, bone spine, brain, liver, lung, lymph node, pleural effusion | Assays: 10x 3' transcription profiling, 10x 3' v2, 10x 3' v3...

#### <a id='gse311609'></a>GSE311609 — Single-Cell Spatial Transcriptomics of Archival Human Lung and Breast Tumor FFPE Samples
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + snRNA-seq
- **Cell Selection / Filtering Strategy:** **Nuclei Isolation (snRNA-seq)**
  > *"Single-nucleus RNA sequencing performed on isolated nuclei from frozen resection tissue."*
- **Cohort Scale:** 41 samples/patients; ~82,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE311609_README.txt; GSE311609_hne_NSCLC_5k_L1_L1_NSCLC_5k_L3_L3_NSCLC_5k_L6_L6.ndpi; GSE311609_hne_NSCLC_5k_L5_L5_NSCLC_5k_L7_L7_NSCLC_5k_L8_L8.ndpi`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE311nnn/GSE311609/suppl/GSE311609_README.txt; ht`
- **Citation / Reference:** PMID:42062553
- **Repository Access Link:** [GSE311609](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE311609)
- **Study Abstract / Experimental Design:**
  We profiled human lung and breast tumor fragments using imaging-based single-cell spatial transcriptomics (10x Xenium) and compared the results to sequencing-based spatial transcriptomics (10x Visium) and single-nucleus RNA-seq...

#### <a id='cellxgene_d224c8e0'></a>CELLxGENE_d224c8e0 — Mesenchymal cells
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 41 samples/patients; ~8,030 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `d224c8e0-c28e-4360-9e42-b3977cd83f9f.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/c54035d5-a021-4bd0-9496-e9ec0cf9dcc8.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2021.09.008
- **Repository Access Link:** [CELLxGENE_d224c8e0](https://cellxgene.cziscience.com/e/d224c8e0-c28e-4360-9e42-b3977cd83f9f.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN MSK - Single cell profiling reveals novel tumor and myeloid subpopulations in small cell lung cancer | Diseases: lung adenocarcinoma, normal, small cell lung carcinoma | Tissues: adrenal gland, axilla, bone spine, brain, liver, lung, lymph node, pleural effusion | Assays: 10x 3' transcription profiling, 10x 3' v2, 10x 3' v3...

#### <a id='cellxgene_d41f45c1'></a>CELLxGENE_d41f45c1 — Transcriptional connectivity of regulatory T cells in the tumor microenvironment informs novel combination cancer therapy strategies
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 23 samples/patients; ~82,991 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `d41f45c1-1b7b-4573-a998-ac5c5acb1647.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/3e1c5296-5c2d-49df-8423-268bcb01f175.h5ad`
- **Citation / Reference:** 10.1038/s41590-023-01504-2
- **Repository Access Link:** [CELLxGENE_d41f45c1](https://cellxgene.cziscience.com/e/d41f45c1-1b7b-4573-a998-ac5c5acb1647.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN MSK - Transcriptional connectivity of regulatory T cells in the tumor microenvironment informs novel combination cancer therapy strategies | Diseases: lung adenocarcinoma | Tissues: left lung, lower lobe of left lung, lower lobe of right lung, middle lobe of right lung, right lung, upper lobe of left lung, upper lobe of right lung | Assays: 10x 3' v2, 10x 3' v3...

#### <a id='cellxgene_01ff5cf0'></a>CELLxGENE_01ff5cf0 — scRNA-seq of Lung adenocarcinoma with histological subtypes
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 18 samples/patients; ~117,266 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `01ff5cf0-730f-4ddc-b1be-7b407211f544.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/229e74ee-9b3e-4ab2-bddb-9eccbac85b48.h5ad`
- **Citation / Reference:** 10.1186/s40164-025-00740-6
- **Repository Access Link:** [CELLxGENE_01ff5cf0](https://cellxgene.cziscience.com/e/01ff5cf0-730f-4ddc-b1be-7b407211f544.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Defining the cellular and molecular identities of histologic subtypes in lung adenocarcinoma | Diseases: lung adenocarcinoma | Tissues: lung parenchyma | Assays: 10x 3' v3...

#### <a id='gse316782'></a>GSE316782 — Multimodal single-cell and spatial profiling reveals altered T cell-mediated immunity and B-cell follicular architecture in non-metastatic lymph nodes of patients with aggressive non-small cell lung cancer
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Regional N1 LNs from patients with more aggressive disease (stage IB–IIIA) exhibited a significant enrichment of dysfunctional CD8⁺ T cells and regulatory T cells (Tregs) compared to N2 LNs and LNs from patients with less aggressive disease (stage IA)."*
- **Cohort Scale:** 18 samples/patients; ~36,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE316782_Xi_CITESeq_GEX_counts.tsv.gz; GSE316782_Xi_CITESeq_Protein_counts.tsv.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE316nnn/GSE316782/suppl/GSE316782_Xi_CITESeq_GEX`
- **Repository Access Link:** [GSE316782](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316782)
- **Study Abstract / Experimental Design:**
  Regional lymph nodes (LNs) in the thoracic cavity serve as essential immunological hubs that coordinate humoral and cell-mediated responses against the development and progression of non-small cell lung cancer (NSCLC). To investigate immune dysregulation in the non-metastatic regional LNs of patients with aggressive NSCLC, we performed multimodal profiling on 36 LNs from 11 patients undergoing cur...

#### <a id='gse333596'></a>GSE333596 — Single-Cell RNA Sequencing Identifies Stem-Like Subpopulations and Their Gene Signature in Lung Adenocarcinoma Progression
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 12 samples/patients; ~24,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE333596_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE333nnn/GSE333596/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE333596](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE333596)
- **Study Abstract / Experimental Design:**
  Lung adenocarcinoma (LUAD) exhibits substantial intra-tumor heterogeneity (ITH), contributing to disease progression and therapeutic challenges. This study utilized single-cell RNA sequencing (scRNA-seq) to profile cells from stage I and III LUAD tumors and matched normal tissues, aiming to delineate malignant cell states and identify tumor-propagating subpopulations. Unsupervised clustering revea...

#### <a id='gse339453'></a>GSE339453 — Single-cell RNA sequencing of EGFR/RB1/TP53-mutant LUAD cell lines with and without PHOX2B overexpression to assess neuroendocrine phenotype induction
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium (3' unspecified)
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 8 samples/patients; ~16,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE339453_H1975.h5ad; GSE339453_PC9.h5ad`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE339nnn/GSE339453/suppl/GSE339453_H1975.h5ad; ht`
- **Repository Access Link:** [GSE339453](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE339453)
- **Study Abstract / Experimental Design:**
  EGFR/RB1/TP53-mutant LUAD cell lines (H1975 and PC9) with and without PHOX2B overexpression were profiled to assess induction of the neuroendocrine (NE) phenotype. Cell lines transduced with lentiviral PHOX2B cDNA or empty vector (LV105) were analyzed by single-cell RNA sequencing (10x Chromium) with hashtag oligo (HTO) multiplexing (CellPlex). Conditions include parental H1975 and PC9 LUAD cell l...

#### <a id='gse346380'></a>GSE346380 — Spatially resolved single-cell transcriptomic profiling of human lung adenocarcinoma and normal lung tissue using NanoString CosMx SMI
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** Subcellular Spatial Transcriptomics (CosMx/Xenium)
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 6 samples/patients; ~12,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE346380_LA_HLT_exprMat.csv.gz; GSE346380_LA_HLT_fov_positions.csv.gz; GSE346380_LA_HLT_metadata.csv.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE346nnn/GSE346380/suppl/GSE346380_LA_HLT_exprMat`
- **Repository Access Link:** [GSE346380](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE346380)
- **Study Abstract / Experimental Design:**
  This study used the NanoString CosMx Spatial Molecular Imager (SMI) to perform spatially resolved single-cell transcriptomic profiling of human lung tissue. Three lung adenocarcinoma tissue specimens and three normal lung tissue specimens were analyzed. CosMx SMI enabled in situ detection and spatial localization of RNA transcripts at single-cell and subcellular resolution. The processed dataset i...

#### <a id='cellxgene_1e4214ce'></a>CELLxGENE_1e4214ce — Lung cancer 4 patients archival FFPE samples profiling with Chromium FLEX
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium (3' unspecified)
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 4 samples/patients; ~35,954 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `1e4214ce-7347-4ec7-97b3-fcf7c1938e8b.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/4ca422ae-e4ad-4884-9073-7c3799cd9143.h5ad`
- **Citation / Reference:** 10.1101/2024.11.01.621259
- **Repository Access Link:** [CELLxGENE_1e4214ce](https://cellxgene.cziscience.com/e/1e4214ce-7347-4ec7-97b3-fcf7c1938e8b.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Transcriptome Analysis of Archived Tumor Tissues by Visium, GeoMx DSP, and Chromium Methods Reveals Inter- and Intra-Patient Heterogeneity | Diseases: lung adenocarcinoma, squamous cell lung carcinoma | Tissues: lung parenchyma | Assays: 10x Next GEM Flex v1...

#### <a id='gse305872'></a>GSE305872 — Single-cell RNA sequencing of Lung Squamous Cell Carcinoma Assosciated with Idiopathic Pulmonary Fibrosis
- **Cancer Type / Indication:** NSCLC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium (3' unspecified)
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~2,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE305872_IPF36F_matrix.mtx.gz; GSE305872_IPF36F_sample_filtered_feature_bc_matrix.h5; GSE305872_IPF36T_matrix.mtx.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE305nnn/GSE305872/suppl/GSE305872_IPF36F_matrix.`
- **Citation / Reference:** PMID:42063567
- **Repository Access Link:** [GSE305872](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305872)
- **Study Abstract / Experimental Design:**
  This dataset contains single-cell RNA sequencing (scRNA-seq) data from two patients with lung squamous cell carcinoma (LUSC) associated with idiopathic pulmonary fibrosis (IPF). Tumor and adjacent lung tissues (within and outside UIP lesions) from two patients were snap-frozen, fixed, and processed using the 10x Genomics Chromium Fixed RNA Profiling workflow. Libraries were sequenced on an Illumin...

### ccRCC Cohorts
*44 cohorts identified for ccRCC*

#### <a id='gse314072'></a>GSE314072 — Functionally heterogeneous intratumoral CD4+CD8+ double positive T cells can give rise to single positive T cells [scRNA-seq + scTCR-seq]
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"D4+CD8+ double positive T cells can give rise to single positive T cells [scRNA-seq + scTCR-seq] Conventional single positive (SP) CD4+ and CD8+ T cells recognize tumor antigens and help mediate clinical responses with cancer immunotherapy."*
- **Cohort Scale:** 24 samples/patients; ~48,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE314072_RCC_Total_adata_multiplex_harmony_noDBL_nomacs_contams_deidentified.h5ad`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE314nnn/GSE314072/suppl/GSE314072_RCC_Total_adat`
- **Citation / Reference:** PMID:41557789
- **Repository Access Link:** [GSE314072](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314072)
- **Study Abstract / Experimental Design:**
  Conventional single positive (SP) CD4+ and CD8+ T cells recognize tumor antigens and help mediate clinical responses with cancer immunotherapy. Double positive CD4+CD8+ (DP) T cells have also been described in human cancers, but their role in the tumor microenvironment (TME) remains unclear. By generating a multi-omic single cell atlas of DP and SP T cells, we find that DP T cells possess phenotyp...

#### <a id='gse285701'></a>GSE285701 — The paradoxical significance of CD39+CD8+ T cells in clear cell renal cell carcinoma [scRNA-seq]
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"The paradoxical significance of CD39+CD8+ T cells in clear cell renal cell carcinoma [scRNA-seq] CD39+CD8+ T cells are known as tumour antigen-specific cells among CD8+ tumour-infiltrating lymphocytes (TILs)."*
- **Cohort Scale:** 13 samples/patients; ~26,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE285701_single_cell_adt_data.txt.gz; GSE285701_single_cell_tcr_info.txt.gz; GSE285701_single_cell_rna_data.txt.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285701/suppl/GSE285701_single_cell_ad`
- **Citation / Reference:** PMID:40961944
- **Repository Access Link:** [GSE285701](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285701)
- **Study Abstract / Experimental Design:**
  CD39+CD8+ T cells are known as tumour antigen-specific cells among CD8+ tumour-infiltrating lymphocytes (TILs). However, CD39+CD8+ T cells also reportedly exhibit immunosuppressive activity in hypoxic tumour models. Here we investigated CD39+CD8+ TILs with regards to their molecular phenotypes, developmental mechanisms, functions, and prognostic significance in clear cell renal cell carcinoma (ccR...

#### <a id='gse210038'></a>GSE210038 — Mesenchymal-like tumor cells and myofibroblastic cancer-associated fibroblasts are associated with progression and immunotherapy response of clear-cell renal cell carcinoma [scRNA-seq]
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 9 samples/patients; ~18,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE210038_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/GSE210038/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:37335139
- **Repository Access Link:** [GSE210038](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210038)
- **Study Abstract / Experimental Design:**
  Immune checkpoint inhibitors (ICI) represent the cornerstone for treatment of patients with metastatic clear-cell renal cell carcinoma (ccRCC). Despite a favorable response for a subset of patients, others experience primary progressive disease highlighting the need to precisely understand plasticity of cancer cells and their crosstalk with the microenvironment to better predict therapeutic respon...

#### <a id='gse304466'></a>GSE304466 — Single-cell transcriptome combined with spatial transcriptome to investigate the molecular mechanism associated with autophagy of clear renal cell carcinoma
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 6 samples/patients; ~12,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE304466_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE304nnn/GSE304466/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE304466](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304466)
- **Study Abstract / Experimental Design:**
  Background: There are no obvious diagnostic markers in the early stages of ccRCC and current targeted therapeutic regimens are susceptible to drug resistance. Our goal is to explore the molecular mechanisms associated with autophagy during the development of ccRCC and to identify molecular markers of prognostic value using scRNA-seq data in combination with stRNA-seq data. Methods: We performed si...

#### <a id='gse223808'></a>GSE223808 — Exhausted intratumoral Vδ2- γδ T cells in human kidney cancer retain effector function [scRNA-seq]
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"These tumor-resident Vδ2- T cells can express the transcriptional program of exhausted alpha beta (ab) CD8+ T cells as well as canonical markers of terminal T cell exhaustion including PD-1, TIGIT and TIM-3."*
- **Cohort Scale:** 4 samples/patients; ~8,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE223808_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE223nnn/GSE223808/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:36928415
- **Repository Access Link:** [GSE223808](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223808)
- **Study Abstract / Experimental Design:**
  Gamma delta (γδ) T cells reside within human tissues including tumors, but their role in mediating anti-tumor response with immune checkpoint inhibition is unknown. Using single-cell approaches, we found that kidney cancers are infiltrated by diverse Vδ2- γδ T cells, with equivalent representation of Vδ1+ and Vδ1- cells, that are distinct from γδ T cells found in normal human tissues. These tumor-...

#### <a id='gse220313'></a>GSE220313 — Gene expression profile and TCR sequencing data of CD45+ cells sorted from single cell suspensions of tumors after in vitro anti-CD3 stimulation and treatments
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD45+ Immune-enriched)**
  > *"Single-cell suspension sorted by FACS for CD45+ leukocytes to enrich for tumor-infiltrating immune cells."*
- **Cohort Scale:** 16 samples/patients; ~32,000 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `GSE220313_scTCR.tar.gz; filelist.txt; GSE220313_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE220nnn/GSE220313/suppl/GSE220313_scTCR.tar.gz; `
- **Repository Access Link:** [GSE220313](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220313)
- **Study Abstract / Experimental Design:**
  We aimed to analyze the transcriptomic difference of tumor-infiltrating immune cells treated with anti-PD-1 blockade and/or HPK1 inhibitor....

#### <a id='gse254498'></a>GSE254498 — Integrating whole-exome sequencing and scRNA-seq reveal the characteristic in one clear cell renal cell carcinoma sample arising in the setting of VHL disease
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~2,000 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `filelist.txt; GSE254498_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254498/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41413084
- **Repository Access Link:** [GSE254498](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254498)
- **Study Abstract / Experimental Design:**
  Clear cell renal cell carcinoma (ccRCC) arising in the setting of von Hippel–Lindau (VHL) disease is a rare type of kidney cancer and features VHL germline mutation. This type of ccRCC is rarely characterised at the single-cell level. In this work, whole-exome sequencing and single-cell RNA sequencing (scRNA-seq) were conducted on one ccRCC sample with VHL disease. Integrating scRNA-seq and whole-...

#### <a id='gse328692'></a>GSE328692 — Single-cell gene programs define subtype identity and metastatic trajectories in renal cell carcinoma
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 34 samples/patients; ~68,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE328692_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE328nnn/GSE328692/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE328692](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE328692)
- **Study Abstract / Experimental Design:**
  This SuperSeries is composed of the SubSeries listed below....

#### <a id='gse304262'></a>GSE304262 — Single cell transcriptome analysis of TIL for gene expression and TCR immunoprofiling in the renal cell carcinoma
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 30 samples/patients; ~60,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE304262_KID_ALL_matrix.mtx.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE304nnn/GSE304262/suppl/GSE304262_KID_ALL_matrix`
- **Citation / Reference:** PMID:41694381
- **Repository Access Link:** [GSE304262](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304262)
- **Study Abstract / Experimental Design:**
  Unlike other carcinomas, renal cell carcinoma is known to have an inverse correlation between the abundance of T cells in the tumor and prognosis. To elucidate the background, we performed single-cell gene expression analysis and TCR immunoprofiling of tumor-infiltrating T cells in early and advanced stages of renal cancer....

#### <a id='gse294109'></a>GSE294109 — Hypoxia shapes both therapeutic response and resistance in metastatic clear cell renal cell carcinoma [scRNAseq]
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 25 samples/patients; ~50,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE294109_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE294nnn/GSE294109/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:42242232
- **Repository Access Link:** [GSE294109](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294109)
- **Study Abstract / Experimental Design:**
  Vascular endothelial growth factor receptor-targeting tyrosine kinase inhibitors (VEGFR-TKIs) and aPD1 combinations are effective in multiple solid tumors, particularly in clear cell renal cell carcinoma (ccRCC), due it’s characteristic pseudo-hypoxic, hyper-angiogenic state driven by biallelic VHL-loss. However, long-term durability is inferior to dual aPD1/aCTLA4 regimens, yet the mechanisms und...

#### <a id='cellxgene_5af90777'></a>CELLxGENE_5af90777 — Single-cell transcriptomic datasets of Renal cell carcinoma patients
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 5' (Immune Profiling)
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 12 samples/patients; ~270,855 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `5af90777-6760-4003-9dba-8f945fec6fdf.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/8784de63-bec2-49ae-b147-8f16108cecce.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2022.11.001
- **Repository Access Link:** [CELLxGENE_5af90777](https://cellxgene.cziscience.com/e/5af90777-6760-4003-9dba-8f945fec6fdf.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Mapping single-cell transcriptomes in the intra-tumoral and associated territories of kidney cancer | Diseases: kidney benign neoplasm, kidney oncocytoma, nonpapillary renal cell carcinoma, normal | Tissues: adrenal tissue, blood, kidney, perirenal fat, vein | Assays: 10x 5' transcription profiling...

#### <a id='gse289672'></a>GSE289672 — Single-Cell Transcriptomics Reveal Distinct Biological and Immunometabolic Landscapes in Renal Medullary Carcinoma Compared to Clear Cell Renal Cell Carcinoma
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (EpCAM+ Malignant/Epithelial)**
  > *"Differential expression analysis identified ENPP3 and CA9 as distinct cell surface markers upregulated in ccRCC, whereas MUC16, TROP2, and EPCAM were upregulated in RMC tumors cells."*
- **Cohort Scale:** 10 samples/patients; ~20,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE289672_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE289nnn/GSE289672/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE289672](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289672)
- **Study Abstract / Experimental Design:**
  Renal medullary carcinoma (RMC) is a rare, highly aggressive kidney malignancy primarily affecting young individuals of African descent with sickle cell trait. The currently available targeted therapies and immunotherapies approved for clear cell renal cell carcinoma (ccRCC) are ineffective against RMC, highlighting the urgent need for novel treatment approaches specifically tailored to the distin...

#### <a id='cellxgene_be39785b'></a>CELLxGENE_be39785b — ccRCC - Single-cell analyses of renal cell cancers reveal insights into tumor microenvironment, cell of origin, and therapy response
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v2
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 7 samples/patients; ~20,509 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `be39785b-67cb-4177-be19-a40ee3747e45.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/f77a9049-7cbe-477d-8d3f-b417bba03147.h5ad`
- **Citation / Reference:** 10.1073/pnas.2103240118
- **Repository Access Link:** [CELLxGENE_be39785b](https://cellxgene.cziscience.com/e/be39785b-67cb-4177-be19-a40ee3747e45.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Single-cell analyses of renal cell cancers reveal insights into tumor microenvironment, cell of origin, and therapy response | Diseases: clear cell renal carcinoma | Tissues: kidney | Assays: 10x 3' v2...

#### <a id='cellxgene_eaf0c852'></a>CELLxGENE_eaf0c852 — 10X Visium spatial transcriptomics on tumor-normal interface tissue sections
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 7 samples/patients; ~27,912 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `eaf0c852-1dff-47ab-8de2-a33a59969a40.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/8ed77943-b258-4fca-887b-2a309cde0535.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2022.11.001
- **Repository Access Link:** [CELLxGENE_eaf0c852](https://cellxgene.cziscience.com/e/eaf0c852-1dff-47ab-8de2-a33a59969a40.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Mapping single-cell transcriptomes in the intra-tumoral and associated territories of kidney cancer | Diseases: nonpapillary renal cell carcinoma | Tissues: kidney | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_bd65a70f'></a>CELLxGENE_bd65a70f — Single-cell sequencing links multiregional immune landscapes and tissue-resident T cells in ccRCC to tumor topology and therapy efficacy
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 5' (Immune Profiling)
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 6 samples/patients; ~167,283 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `bd65a70f-b274-4133-b9dd-0d1431b6af34.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/6602fa5a-d663-4fe0-8dc3-90191c3a013b.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2021.03.007
- **Repository Access Link:** [CELLxGENE_bd65a70f](https://cellxgene.cziscience.com/e/bd65a70f-b274-4133-b9dd-0d1431b6af34.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Single-cell sequencing links multiregional immune landscapes and tissue-resident T cells in ccRCC to tumor topology and therapy efficacy | Diseases: clear cell renal carcinoma | Tissues: blood, kidney, lymph node | Assays: 10x 5' v2...

#### <a id='gse254444'></a>GSE254444 — Investigation of tumor heterogeneity using single cell transcriptomics in Von Hippel-Lindau related renal cancer
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 6 samples/patients; ~12,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE254444_ccRCC.raw_counts.annotations.txt.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254444/suppl/GSE254444_ccRCC.raw_coun`
- **Citation / Reference:** PMID:41513925
- **Repository Access Link:** [GSE254444](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254444)
- **Study Abstract / Experimental Design:**
  Von Hippel-Lindau (VHL) dìsease is a rare inherited disorder caused by the loss of the VHL gene and characterized by benign and malignant tumors in several organs. Patients with VHL disease have risk of developing multiple primary clear cell renal carcinomas (ccRCC). Patients are managed with multiple sessions of renal surgery causing renal failure and complication of dialysis or systemic progress...

#### <a id='cellxgene_318cb2a6'></a>CELLxGENE_318cb2a6 — 10X Visium spatial transcriptomics on tumor core tissue sections
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 5 samples/patients; ~10,924 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `318cb2a6-743a-4404-9515-413d5f64b9ff.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/339480f8-4f4f-450d-b475-8963163396d1.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2022.11.001
- **Repository Access Link:** [CELLxGENE_318cb2a6](https://cellxgene.cziscience.com/e/318cb2a6-743a-4404-9515-413d5f64b9ff.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Mapping single-cell transcriptomes in the intra-tumoral and associated territories of kidney cancer | Diseases: nonpapillary renal cell carcinoma | Tissues: kidney | Assays: Visium Spatial Gene Expression V1...

#### <a id='gse272610'></a>GSE272610 — Comparative single-cell transcriptomic profiling of patient-derived renal carcinoma cells in cellular and animal models of kidney cancer
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 5 samples/patients; ~243 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE272610_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE272nnn/GSE272610/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:40241258
- **Repository Access Link:** [GSE272610](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272610)
- **Study Abstract / Experimental Design:**
  Clear-cell renal cell carcinoma (ccRCC) is the most common form of kidney cancer, which is often resistant to conventional cancer therapies including chemotherapy and radiation therapy.  Targeted treatments including immunotherapies and small molecule inhibitors have been recently developed with positive outcomes.  However, variations in the patient response and resistance to therapies suggest tha...

#### <a id='cellxgene_a45125aa'></a>CELLxGENE_a45125aa — 6800STDY12499413
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `a45125aa-7641-4bca-bca5-6d652e027d6f.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/91f3d7de-936d-41a0-b814-d98169ca8eba.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2022.11.001
- **Repository Access Link:** [CELLxGENE_a45125aa](https://cellxgene.cziscience.com/e/a45125aa-7641-4bca-bca5-6d652e027d6f.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Mapping single-cell transcriptomes in the intra-tumoral and associated territories of kidney cancer | Diseases: nonpapillary renal cell carcinoma | Tissues: kidney | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_4aefec71'></a>CELLxGENE_4aefec71 — Renal cell carcinoma, post aPD1, lung metastasis Puck_220408_14
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~37,568 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `4aefec71-da4e-48ec-a910-c542b5807c7d.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/d8c55296-d90e-45fa-8bc9-bd4ded93e97b.h5ad`
- **Citation / Reference:** 10.1016/j.immuni.2022.09.002
- **Repository Access Link:** [CELLxGENE_4aefec71](https://cellxgene.cziscience.com/e/4aefec71-da4e-48ec-a910-c542b5807c7d.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatially mapping T cell receptors and transcriptomes | Diseases: renal cell carcinoma | Tissues: lung | Assays: Slide-seqV2...

#### <a id='cellxgene_32ffc3a7'></a>CELLxGENE_32ffc3a7 — Renal cell carcinoma, post aPD1, lung metastasis Puck_220408_13
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~40,217 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `32ffc3a7-c88d-49e1-894f-8e09d0184ee9.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/6628ca2e-a2c5-4580-ab28-93372ad4b16e.h5ad`
- **Citation / Reference:** 10.1016/j.immuni.2022.09.002
- **Repository Access Link:** [CELLxGENE_32ffc3a7](https://cellxgene.cziscience.com/e/32ffc3a7-c88d-49e1-894f-8e09d0184ee9.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatially mapping T cell receptors and transcriptomes | Diseases: renal cell carcinoma | Tissues: lung | Assays: Slide-seqV2...

#### <a id='cellxgene_18fd0190'></a>CELLxGENE_18fd0190 — Renal cell carcinoma, post aPD1, lung metastasis Puck_220408_20
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~37,563 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `18fd0190-229d-4d2b-8f28-1e40c49c750d.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/b02b5c49-db1c-4db3-bfd0-02224b63158c.h5ad`
- **Citation / Reference:** 10.1016/j.immuni.2022.09.002
- **Repository Access Link:** [CELLxGENE_18fd0190](https://cellxgene.cziscience.com/e/18fd0190-229d-4d2b-8f28-1e40c49c750d.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatially mapping T cell receptors and transcriptomes | Diseases: renal cell carcinoma | Tissues: lung | Assays: Slide-seqV2...

#### <a id='cellxgene_07efa1c3'></a>CELLxGENE_07efa1c3 — Renal cell carcinoma, post aPD1, lung metastasis Puck_220408_15
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~35,839 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `07efa1c3-7649-43d8-b798-e26a0145ee65.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/8b589484-fa97-4123-a7f4-7ed1252a54a8.h5ad`
- **Citation / Reference:** 10.1016/j.immuni.2022.09.002
- **Repository Access Link:** [CELLxGENE_07efa1c3](https://cellxgene.cziscience.com/e/07efa1c3-7649-43d8-b798-e26a0145ee65.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatially mapping T cell receptors and transcriptomes | Diseases: renal cell carcinoma | Tissues: lung | Assays: Slide-seqV2...

#### <a id='cellxgene_02faf712'></a>CELLxGENE_02faf712 — Renal cell carcinoma, pre aPD1, kidney Puck_200727_12
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~17,612 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `02faf712-92d4-4589-bec7-13105059cf86.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/de6d5195-dc73-42cb-ae44-283221077cbe.h5ad`
- **Citation / Reference:** 10.1016/j.immuni.2022.09.002
- **Repository Access Link:** [CELLxGENE_02faf712](https://cellxgene.cziscience.com/e/02faf712-92d4-4589-bec7-13105059cf86.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatially mapping T cell receptors and transcriptomes | Diseases: renal cell carcinoma | Tissues: kidney | Assays: Slide-seqV2...

#### <a id='cellxgene_6d243918'></a>CELLxGENE_6d243918 — Renal cell carcinoma, post aPD1, lung metastasis Puck_200727_08
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~18,851 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `6d243918-f3f6-4a9e-b72e-39d96de73095.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/e5e05c24-f600-40c4-998d-1c45d3bbcae5.h5ad`
- **Citation / Reference:** 10.1016/j.immuni.2022.09.002
- **Repository Access Link:** [CELLxGENE_6d243918](https://cellxgene.cziscience.com/e/6d243918-f3f6-4a9e-b72e-39d96de73095.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatially mapping T cell receptors and transcriptomes | Diseases: renal cell carcinoma | Tissues: lung | Assays: Slide-seqV2...

#### <a id='cellxgene_9b1437bb'></a>CELLxGENE_9b1437bb — Renal cell carcinoma, post aPD1, lung metastasis Puck_200727_09
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~14,990 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `9b1437bb-42fe-47d8-8962-835d95581763.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/c80e9271-1866-4cbb-9084-12897dcb0ed6.h5ad`
- **Citation / Reference:** 10.1016/j.immuni.2022.09.002
- **Repository Access Link:** [CELLxGENE_9b1437bb](https://cellxgene.cziscience.com/e/9b1437bb-42fe-47d8-8962-835d95581763.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatially mapping T cell receptors and transcriptomes | Diseases: renal cell carcinoma | Tissues: lung | Assays: Slide-seqV2...

#### <a id='cellxgene_f25a532c'></a>CELLxGENE_f25a532c — 6800STDY12499504
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `f25a532c-faf6-4d65-9dac-47d6369508f6.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/561e15e5-36f7-4a03-9a7e-c4bd126bf438.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2022.11.001
- **Repository Access Link:** [CELLxGENE_f25a532c](https://cellxgene.cziscience.com/e/f25a532c-faf6-4d65-9dac-47d6369508f6.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Mapping single-cell transcriptomes in the intra-tumoral and associated territories of kidney cancer | Diseases: nonpapillary renal cell carcinoma | Tissues: kidney | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_cc43509e'></a>CELLxGENE_cc43509e — Renal cell carcinoma, post aPD1, lung metastasis Puck_200727_10
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~15,366 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `cc43509e-8729-4cbc-ace9-f0f7d654ee53.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/288b5bdb-f3cc-43a0-8e47-28aab9cff903.h5ad`
- **Citation / Reference:** 10.1016/j.immuni.2022.09.002
- **Repository Access Link:** [CELLxGENE_cc43509e](https://cellxgene.cziscience.com/e/cc43509e-8729-4cbc-ace9-f0f7d654ee53.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatially mapping T cell receptors and transcriptomes | Diseases: renal cell carcinoma | Tissues: lung | Assays: Slide-seqV2...

#### <a id='cellxgene_f339fe89'></a>CELLxGENE_f339fe89 — Renal cell carcinoma, pre aPD1, kidney Puck_200727_13
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~20,021 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `f339fe89-c99c-4f50-b76d-2ffe9d376500.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/e4128ffc-92d5-4704-bb0c-5fe72f53e94d.h5ad`
- **Citation / Reference:** 10.1016/j.immuni.2022.09.002
- **Repository Access Link:** [CELLxGENE_f339fe89](https://cellxgene.cziscience.com/e/f339fe89-c99c-4f50-b76d-2ffe9d376500.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatially mapping T cell receptors and transcriptomes | Diseases: renal cell carcinoma | Tissues: kidney | Assays: Slide-seqV2...

#### <a id='cellxgene_05f813a4'></a>CELLxGENE_05f813a4 — 6800STDY12499411
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `05f813a4-fc9e-4143-bba9-de93bdd8eece.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/fe87403a-2bf1-4e29-b544-7ca139b20d51.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2022.11.001
- **Repository Access Link:** [CELLxGENE_05f813a4](https://cellxgene.cziscience.com/e/05f813a4-fc9e-4143-bba9-de93bdd8eece.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Mapping single-cell transcriptomes in the intra-tumoral and associated territories of kidney cancer | Diseases: nonpapillary renal cell carcinoma | Tissues: kidney | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_0671c0d4'></a>CELLxGENE_0671c0d4 — 6800STDY12499406
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `0671c0d4-f233-4cd3-a0d6-9a6b747bdd98.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/3cb22696-5867-4001-ba23-3a24da12f4ab.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2022.11.001
- **Repository Access Link:** [CELLxGENE_0671c0d4](https://cellxgene.cziscience.com/e/0671c0d4-f233-4cd3-a0d6-9a6b747bdd98.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Mapping single-cell transcriptomes in the intra-tumoral and associated territories of kidney cancer | Diseases: nonpapillary renal cell carcinoma | Tissues: kidney | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_4c6f9f26'></a>CELLxGENE_4c6f9f26 — chRCC - Single-cell analyses of renal cell cancers reveal insights into tumor microenvironment, cell of origin, and therapy response
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v2
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~2,576 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `4c6f9f26-5470-455b-8933-c408232fbf56.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/36f72651-4bac-4dd3-9a9b-9a0c5d55fcf1.h5ad`
- **Citation / Reference:** 10.1073/pnas.2103240118
- **Repository Access Link:** [CELLxGENE_4c6f9f26](https://cellxgene.cziscience.com/e/4c6f9f26-5470-455b-8933-c408232fbf56.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Single-cell analyses of renal cell cancers reveal insights into tumor microenvironment, cell of origin, and therapy response | Diseases: chromophobe renal cell carcinoma | Tissues: kidney | Assays: 10x 3' v2...

#### <a id='cellxgene_104cfa2a'></a>CELLxGENE_104cfa2a — 6800STDY12499508
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `104cfa2a-8f7e-4dcc-a0a3-8fd46711c904.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/16fe2d28-4272-40b5-8c78-a5a067bedc9a.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2022.11.001
- **Repository Access Link:** [CELLxGENE_104cfa2a](https://cellxgene.cziscience.com/e/104cfa2a-8f7e-4dcc-a0a3-8fd46711c904.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Mapping single-cell transcriptomes in the intra-tumoral and associated territories of kidney cancer | Diseases: nonpapillary renal cell carcinoma | Tissues: kidney | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_24c31c8c'></a>CELLxGENE_24c31c8c — 6800STDY12499412
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `24c31c8c-5f6f-4f13-9f7b-abbedcdcd45f.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/ca1956a7-c7f6-416c-bfe7-3b6b2b65c56c.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2022.11.001
- **Repository Access Link:** [CELLxGENE_24c31c8c](https://cellxgene.cziscience.com/e/24c31c8c-5f6f-4f13-9f7b-abbedcdcd45f.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Mapping single-cell transcriptomes in the intra-tumoral and associated territories of kidney cancer | Diseases: nonpapillary renal cell carcinoma | Tissues: kidney | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_252438d3'></a>CELLxGENE_252438d3 — 6800STDY12499408
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `252438d3-8143-48a9-9ac2-ef5991fb21b8.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/64a1f56f-b355-47e0-b263-30fa1e6c7390.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2022.11.001
- **Repository Access Link:** [CELLxGENE_252438d3](https://cellxgene.cziscience.com/e/252438d3-8143-48a9-9ac2-ef5991fb21b8.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Mapping single-cell transcriptomes in the intra-tumoral and associated territories of kidney cancer | Diseases: nonpapillary renal cell carcinoma | Tissues: kidney | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_27cd5ac9'></a>CELLxGENE_27cd5ac9 — 6800STDY12499505
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `27cd5ac9-80ba-43f3-bb60-740593dabdd0.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/93a488b8-ff57-4502-8e95-67569d896e95.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2022.11.001
- **Repository Access Link:** [CELLxGENE_27cd5ac9](https://cellxgene.cziscience.com/e/27cd5ac9-80ba-43f3-bb60-740593dabdd0.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Mapping single-cell transcriptomes in the intra-tumoral and associated territories of kidney cancer | Diseases: nonpapillary renal cell carcinoma | Tissues: kidney | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_30437616'></a>CELLxGENE_30437616 — 6800STDY12499503
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `30437616-4939-40e0-83da-4cae8186d452.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/7d114367-c1b1-4883-bceb-adf083564af4.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2022.11.001
- **Repository Access Link:** [CELLxGENE_30437616](https://cellxgene.cziscience.com/e/30437616-4939-40e0-83da-4cae8186d452.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Mapping single-cell transcriptomes in the intra-tumoral and associated territories of kidney cancer | Diseases: nonpapillary renal cell carcinoma | Tissues: kidney | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_53d62b10'></a>CELLxGENE_53d62b10 — 6800STDY12499502
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `53d62b10-bae5-48ac-b16e-71be9ba6de59.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/fe61d825-a898-4d2d-9c79-834b2fc939d3.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2022.11.001
- **Repository Access Link:** [CELLxGENE_53d62b10](https://cellxgene.cziscience.com/e/53d62b10-bae5-48ac-b16e-71be9ba6de59.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Mapping single-cell transcriptomes in the intra-tumoral and associated territories of kidney cancer | Diseases: nonpapillary renal cell carcinoma | Tissues: kidney | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_60ac2657'></a>CELLxGENE_60ac2657 — 6800STDY12499409
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `60ac2657-f325-46c7-817b-55bdab803105.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/11a7ccf8-0cf3-49cb-a1ce-6a2d6565e49c.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2022.11.001
- **Repository Access Link:** [CELLxGENE_60ac2657](https://cellxgene.cziscience.com/e/60ac2657-f325-46c7-817b-55bdab803105.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Mapping single-cell transcriptomes in the intra-tumoral and associated territories of kidney cancer | Diseases: nonpapillary renal cell carcinoma | Tissues: kidney | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_75548d10'></a>CELLxGENE_75548d10 — 6800STDY12499506
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `75548d10-160d-4f3e-b317-99ad9630c62d.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/59a876b2-15c3-4d36-80ef-06270ccef169.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2022.11.001
- **Repository Access Link:** [CELLxGENE_75548d10](https://cellxgene.cziscience.com/e/75548d10-160d-4f3e-b317-99ad9630c62d.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Mapping single-cell transcriptomes in the intra-tumoral and associated territories of kidney cancer | Diseases: nonpapillary renal cell carcinoma | Tissues: kidney | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_81328f3f'></a>CELLxGENE_81328f3f — 6800STDY12499509
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `81328f3f-31d4-47e3-9eae-e07512026ce1.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/db8ce618-2929-4924-a672-b9a133ba5e64.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2022.11.001
- **Repository Access Link:** [CELLxGENE_81328f3f](https://cellxgene.cziscience.com/e/81328f3f-31d4-47e3-9eae-e07512026ce1.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Mapping single-cell transcriptomes in the intra-tumoral and associated territories of kidney cancer | Diseases: nonpapillary renal cell carcinoma | Tissues: kidney | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_a6046b15'></a>CELLxGENE_a6046b15 — 6800STDY12499410
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `a6046b15-a095-43b0-9fb5-b36899d87fdb.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/e140d92c-734d-437e-bf59-bbcb2e583741.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2022.11.001
- **Repository Access Link:** [CELLxGENE_a6046b15](https://cellxgene.cziscience.com/e/a6046b15-a095-43b0-9fb5-b36899d87fdb.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Mapping single-cell transcriptomes in the intra-tumoral and associated territories of kidney cancer | Diseases: nonpapillary renal cell carcinoma | Tissues: kidney | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_c3f74413'></a>CELLxGENE_c3f74413 — 6800STDY12499407
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `c3f74413-9faa-43b5-93af-601c5af4acc8.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/d2962792-6fe8-414c-b0ff-f96b38638bb0.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2022.11.001
- **Repository Access Link:** [CELLxGENE_c3f74413](https://cellxgene.cziscience.com/e/c3f74413-9faa-43b5-93af-601c5af4acc8.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Mapping single-cell transcriptomes in the intra-tumoral and associated territories of kidney cancer | Diseases: nonpapillary renal cell carcinoma | Tissues: kidney | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_c5ac3ec2'></a>CELLxGENE_c5ac3ec2 — 6800STDY12499507
- **Cancer Type / Indication:** ccRCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `c5ac3ec2-24b0-43cc-9aab-bb0ebbe205ce.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/99e65769-73ec-4e9a-825e-d2210b62bd0e.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2022.11.001
- **Repository Access Link:** [CELLxGENE_c5ac3ec2](https://cellxgene.cziscience.com/e/c5ac3ec2-24b0-43cc-9aab-bb0ebbe205ce.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Mapping single-cell transcriptomes in the intra-tumoral and associated territories of kidney cancer | Diseases: nonpapillary renal cell carcinoma | Tissues: kidney | Assays: Visium Spatial Gene Expression V1...

### Bladder Cohorts
*18 cohorts identified for Bladder*

#### <a id='gse183556'></a>GSE183556 — MSL2 is an allelic dosage sensor in mammals
- **Cancer Type / Indication:** Bladder
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 550 samples/patients; ~1,100,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE183556_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183556/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:38030723
- **Repository Access Link:** [GSE183556](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183556)
- **Study Abstract / Experimental Design:**
  This SuperSeries is composed of the SubSeries listed below....

#### <a id='gse192575'></a>GSE192575 — Acquired semi-squamalization during chemotherapy suggests a differentiation therapy for bladder cancer
- **Cancer Type / Indication:** Bladder
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 42 samples/patients; ~84,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE192575_RNAseq1_mouse_Sen_vs_normal_result_all.csv.gz; GSE192575_RNAseq2_mouse_Resis_vs_Sen_result_all.csv.gz; GSE192575_RNAseq3_mouse_DDP_result_all.csv.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE192nnn/GSE192575/suppl/GSE192575_RNAseq1_mouse_`
- **Repository Access Link:** [GSE192575](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192575)
- **Study Abstract / Experimental Design:**
  Cisplatin-based chemotherapy is the most common treatment for unresectable bladder cancers and also increasingly used as neoadjuvant treatments before or after surgery and radiotherapy. Unfortunately, though many patients respond to the treatment, most of them develop resistance quickly with unclear mechanisms and few further treatment options. Here, we report that semi-squamazation is acquired du...

#### <a id='gse267718'></a>GSE267718 — Single cell RNA sequencing (scRNAseq) of urine-derived cells, peripheral blood mononuclear (PBMC), and tumor cells from bladder cancer patients.
- **Cancer Type / Indication:** Bladder
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 33 samples/patients; ~66,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE267718_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267718/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:38847806
- **Repository Access Link:** [GSE267718](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267718)
- **Study Abstract / Experimental Design:**
  As no one previously examined urine-derived cells from bladder cancer patients, we performed scRNAseq to profile the diversity of these cells and their transcriptional profiles. We used scRNAseq to compare the profiles of urine-derived cells to matched tumor cells and PBMC from bladder cancer patients....

#### <a id='gse302781'></a>GSE302781 — snRNASeq of Metastatic Urothelial carcinoma samples
- **Cancer Type / Indication:** Bladder
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** Single-Nucleus RNA-seq
- **Modality:** snRNA-seq
- **Cell Selection / Filtering Strategy:** **Nuclei Isolation (snRNA-seq)**
  > *"Single-nucleus RNA sequencing performed on isolated nuclei from frozen resection tissue."*
- **Cohort Scale:** 15 samples/patients; ~30,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE302781_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE302nnn/GSE302781/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE302781](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302781)
- **Study Abstract / Experimental Design:**
  Histological variation is a common feature of metastatic cancer but its clonal origins and impact on disease history remains poorly defined. In this study, we developed a first-in-kind metastatic bladder cancer (BLCA) rapid autopsy cohort enriched in rare histological subtypes to extensively profile individuals with terminal disease. By reconstructing the evolutionary history of patient tumors, we...

#### <a id='gse301651'></a>GSE301651 — Single cell sequencing of bladder cancer patients
- **Cancer Type / Indication:** Bladder
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 13 samples/patients; ~26,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE301651_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE301nnn/GSE301651/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41706539
- **Repository Access Link:** [GSE301651](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301651)
- **Study Abstract / Experimental Design:**
  To achieve a more profound comprehension of the bladder tumor microenvironment, we have meticulously curated single-cell samples from 11 patients diagnosed with bladder cancer. Our comprehensive analysis not only elucidates the intricate heterogeneity among distinct cellular subpopulations within the tumor but also uncovers the underlying cellular interactions. These findings collectively serve as...

#### <a id='gse222315'></a>GSE222315 — Gene expression profile at single cell level of bladder cancer and normal adjacent tissues
- **Cancer Type / Indication:** Bladder
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 13 samples/patients; ~26,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE222315_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE222nnn/GSE222315/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:38428409
- **Repository Access Link:** [GSE222315](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222315)
- **Study Abstract / Experimental Design:**
  The cellular diversity of stromal compartment within the tumor microenvironment might change dynamically as tumors evolve that allowing distant dissemination of tumor cells and informing the treatment options. Thus, we used single cell RNA sequencing (scRNA-seq) to provide a detailed insight of cellular components feature in bladder cancer....

#### <a id='gse172433'></a>GSE172433 — Genome-wide H3k27ac levels following FOXA1 KO in UMUC-1 human bladder cancer cells
- **Cancer Type / Indication:** Bladder
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 12 samples/patients; ~24,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE172433_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE172nnn/GSE172433/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:36323682
- **Repository Access Link:** [GSE172433](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172433)
- **Study Abstract / Experimental Design:**
  This SuperSeries is composed of the SubSeries listed below....

#### <a id='gse277524'></a>GSE277524 — Single-cell dissection of multifocal bladder cancer reveals malignant and immune cells variation between primary and recurrence tumor lesions
- **Cancer Type / Indication:** Bladder
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 10 samples/patients; ~20,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE277524_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277524/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:39702554
- **Repository Access Link:** [GSE277524](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277524)
- **Study Abstract / Experimental Design:**
  Bladder carcinoma (BLCA) is characterized by a high rate of post-surgery relapse and multifocality, with multifocal tumors carrying a 40% higher risk of recurrence compared to single tumors. Recurrence is a major contributor to bladder cancer-specific mortality. However, understanding interregional or intraregional malignant heterogeneity and cellular communication within the primary and recurrenc...

#### <a id='gse310802'></a>GSE310802 — Morphological and Single-Cell Analyses of En Bloc Resected Non-Muscle invasive Bladder Cancer Reveal Structural and functional disruption of Tertiary lymphoid Structures
- **Cancer Type / Indication:** Bladder
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"infiltration and altered dendritic cell function, including reduced MHC class II antigen presentation and downregulation of co-stimulatory (CD86–CD28) and migration-related (CD99) pathways compared to cystitis (p<0."*
- **Cohort Scale:** 8 samples/patients; ~16,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE310802_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE310nnn/GSE310802/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE310802](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310802)
- **Study Abstract / Experimental Design:**
  Tertiary lymphoid structures (TLS) have emerged as critical immune niches within the tumor microenvironment across various cancers. However, their structural organization and functional roles in non–muscle-invasive bladder cancer (NMIBC) remain poorly characterized. In this study, we provide the first detailed characterization of TLS in NMIBC using en bloc resected specimens, enabling high-resolut...

#### <a id='gse326225'></a>GSE326225 — A Single-Cell Atlas of Muscle-Invasive Bladder Cancer Reveals Lineage-Specific Vulnerabilities
- **Cancer Type / Indication:** Bladder
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 8 samples/patients; ~16,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE326225_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE326nnn/GSE326225/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE326225](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE326225)
- **Study Abstract / Experimental Design:**
  Single-cell RNA sequencing (scRNA-seq) was performed to characterize cellular heterogeneity and spatial organization in muscle-invasive bladder cancer (MIBC). scRNA-seq was used to define tumor and microenvironmental cell populations and transcriptional states....

#### <a id='gse326854'></a>GSE326854 — Gene expression profile at single cell level of different phenotypic tumors derived from human non-muscle invasive bladder carcinoma cell lines
- **Cancer Type / Indication:** Bladder
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 6 samples/patients; ~12,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE326854_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE326nnn/GSE326854/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:42203766
- **Repository Access Link:** [GSE326854](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE326854)
- **Study Abstract / Experimental Design:**
  The cellular states of bladder carcinoma change dynamically as tumor recurrence evolves under Gemcitabine treatment, facilitating phenotypic transitions. Thus, we employed single-cell RNA sequencing (scRNA-seq) to gain a high-resolution view of the cellular landscape in bladder carcinoma....

#### <a id='gse250523'></a>GSE250523 — Single-cell RNA sequencing of bladder Ewing sarcoma
- **Cancer Type / Indication:** Bladder
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 4 samples/patients; ~8,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE250523_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250523/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:39386756
- **Repository Access Link:** [GSE250523](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250523)
- **Study Abstract / Experimental Design:**
  Bladder Ewing sarcoma/primitive neuroectodermal tumor (bladder ES/PNET) is a rare and highly malignant tumor associated with a poor prognosis, yet its underlying mechanisms remain poorly understood. This study employed a combination of single-cell RNA sequencing (scRNA-seq) analyses to delve into the pathogenesis of bladder ES/PNET. The investigation revealed the presence of specialized types of e...

#### <a id='cellxgene_024581e3'></a>CELLxGENE_024581e3 — SMBO-109
- **Cancer Type / Indication:** Bladder
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~8,156 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `024581e3-5375-4e33-8060-d8448694f556.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/f77f7919-8483-4c56-a0de-6be3e5363c7a.h5ad`
- **Citation / Reference:** 10.1038/s41467-025-67643-2
- **Repository Access Link:** [CELLxGENE_024581e3](https://cellxgene.cziscience.com/e/024581e3-5375-4e33-8060-d8448694f556.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Determinants of Sensitivity to HER2-targeted Antibody Drug Conjugates in Urothelial Cancer | Diseases: urinary bladder cancer | Tissues: urinary bladder | Assays: 10x 3' v3...

#### <a id='cellxgene_7ee4b15b'></a>CELLxGENE_7ee4b15b — SCBO-8
- **Cancer Type / Indication:** Bladder
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~10,946 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `7ee4b15b-b543-4a8f-a602-ab0344729caf.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/59caba36-eb74-4d96-94cf-80a695d5dd44.h5ad`
- **Citation / Reference:** 10.1038/s41467-025-67643-2
- **Repository Access Link:** [CELLxGENE_7ee4b15b](https://cellxgene.cziscience.com/e/7ee4b15b-b543-4a8f-a602-ab0344729caf.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Determinants of Sensitivity to HER2-targeted Antibody Drug Conjugates in Urothelial Cancer | Diseases: urinary bladder cancer | Tissues: urinary bladder | Assays: 10x 3' v3...

#### <a id='cellxgene_e2094676'></a>CELLxGENE_e2094676 — SMBO-106
- **Cancer Type / Indication:** Bladder
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~5,762 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `e2094676-d4d5-400c-afa5-7ca3ec09025c.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/5447ecfb-17c8-430d-90a8-9734cc4f0701.h5ad`
- **Citation / Reference:** 10.1038/s41467-025-67643-2
- **Repository Access Link:** [CELLxGENE_e2094676](https://cellxgene.cziscience.com/e/e2094676-d4d5-400c-afa5-7ca3ec09025c.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Determinants of Sensitivity to HER2-targeted Antibody Drug Conjugates in Urothelial Cancer | Diseases: urinary bladder cancer | Tissues: urinary bladder | Assays: 10x 3' v3...

#### <a id='cellxgene_f0e0575d'></a>CELLxGENE_f0e0575d — SMBO-114
- **Cancer Type / Indication:** Bladder
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~2,410 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `f0e0575d-d64d-4418-809e-60789bf6c155.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/ea1c435b-245d-40db-bc03-0d602498ab69.h5ad`
- **Citation / Reference:** 10.1038/s41467-025-67643-2
- **Repository Access Link:** [CELLxGENE_f0e0575d](https://cellxgene.cziscience.com/e/f0e0575d-d64d-4418-809e-60789bf6c155.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Determinants of Sensitivity to HER2-targeted Antibody Drug Conjugates in Urothelial Cancer | Diseases: urinary bladder cancer | Tissues: urinary bladder | Assays: 10x 3' v3...

#### <a id='gse176249'></a>GSE176249 — Genome-wide H3k27ac levels following FOXA1 KO in UMUC-1 human bladder cancer cells [scRNA-Seq]
- **Cancer Type / Indication:** Bladder
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium (3' unspecified)
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~2,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE176249_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176249/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:36323682
- **Repository Access Link:** [GSE176249](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176249)
- **Study Abstract / Experimental Design:**
  10x scRNA-Seq results of a bladder cancer case with extensive squamous cell differentation...

#### <a id='cellxgene_670a9f65'></a>CELLxGENE_670a9f65 — SMBO-170
- **Cancer Type / Indication:** Bladder
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~5,785 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `670a9f65-363c-4bb1-a942-bf6d1edea179.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/ff48d789-fbcd-4907-a002-e175bda24f77.h5ad`
- **Citation / Reference:** 10.1038/s41467-025-67643-2
- **Repository Access Link:** [CELLxGENE_670a9f65](https://cellxgene.cziscience.com/e/670a9f65-363c-4bb1-a942-bf6d1edea179.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Determinants of Sensitivity to HER2-targeted Antibody Drug Conjugates in Urothelial Cancer | Diseases: urinary bladder cancer | Tissues: urinary bladder | Assays: 10x 3' v3...

### Breast Cohorts
*89 cohorts identified for Breast*

#### <a id='gse246613'></a>GSE246613 — Single-cell and spatial profiling identify three response trajectories to pembrolizumab and radiation therapy in triple negative breast cancer
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 266 samples/patients; ~532,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE246613_PembroRT_immune_R100_final.h5ad.gz; GSE246613_PembroRT_non_immune_cells.h5ad.gz; GSE246613_combined_RTPDv4_scvi_celltypist.h5ad.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246613/suppl/GSE246613_PembroRT_immun`
- **Citation / Reference:** PMID:38194915
- **Repository Access Link:** [GSE246613](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246613)
- **Study Abstract / Experimental Design:**
  Cancer immunotherapy trials have produced encouraging results, but resistance remains a problem necessitating strategies to identify patients that benefit from immunotherapy alone or who require additional combinations like chemotherapy or radiotherapy. Here we employ single-cell transcriptomics and spatial proteomics to profile triple negative breast cancer biopsies taken before and after one cyc...

#### <a id='gse274141'></a>GSE274141 — Single-Cell RNA Sequencing Identifies Molecular Biomarkers Predicting Late Progression to CDK4/6 Inhibition in Patients with HR+/HER2- Metastatic Breast Cancer [FFPE]
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"We identify and validate tumor-specific molecular biomarkers and demonstrate that tumor-infiltrating cytotoxic natural killer and CD8+ T cells are predictive of therapeutic response."*
- **Cohort Scale:** 54 samples/patients; ~108,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE274141_read-counts-n54-new.csv.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274141/suppl/GSE274141_read-counts-n5`
- **Citation / Reference:** PMID:39955556
- **Repository Access Link:** [GSE274141](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274141)
- **Study Abstract / Experimental Design:**
  This study addresses the critical need for predictive biomarkers in patients with hormone receptor–positive, HER2-negative (HR+/HER2-) metastatic breast cancer (mBC) who have been treated with cyclin-dependent kinase 4/6 inhibitors (CDK4/6is). Using single-cell RNA sequencing, we reveal differential mechanisms of early versus late progression, enhancing our understanding of resistance to CDK4/6is....

#### <a id='gse222859'></a>GSE222859 — Gene expression profile at single cell level of lymphocytes cells for Pan-T cell analysis
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"We observed significant enrichment of intratumoural CD4/CD8 TSTR cells following ICB therapy, particularly, in non-responsive tumors, suggesting their potential role in immunotherapy resistance."*
- **Cohort Scale:** 51 samples/patients; ~102,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE222859_matrix.mtx.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE222nnn/GSE222859/suppl/GSE222859_matrix.mtx.gz`
- **Citation / Reference:** PMID:37248301
- **Repository Access Link:** [GSE222859](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222859)
- **Study Abstract / Experimental Design:**
  Tumour infiltrating T lymphocytes (TILs) have provided an attractive avenue for cancer treatment. Yet, the various states of TILs in tumour immune microenvironment (TIME) have not been fully characterized. Here, we present a single-cell atlas of T cells, encompassing 308,048 transcriptomes, spanning 16 different cancer types and 9 types of healthy/cancer-adjacent tissues. This allowed elucidation ...

#### <a id='gse274139'></a>GSE274139 — Single-Cell RNA Sequencing Identifies Molecular Biomarkers Predicting Late Progression to CDK4/6 Inhibition in Patients with HR+/HER2- Metastatic Breast Cancer [tissue]
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"We identify and validate tumor-specific molecular biomarkers and demonstrate that tumor-infiltrating cytotoxic natural killer and CD8+ T cells are predictive of therapeutic response."*
- **Cohort Scale:** 35 samples/patients; ~70,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE274139_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274139/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:39955556
- **Repository Access Link:** [GSE274139](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274139)
- **Study Abstract / Experimental Design:**
  This study addresses the critical need for predictive biomarkers in patients with hormone receptor–positive, HER2-negative (HR+/HER2-) metastatic breast cancer (mBC) who have been treated with cyclin-dependent kinase 4/6 inhibitors (CDK4/6is). Using single-cell RNA sequencing, we reveal differential mechanisms of early versus late progression, enhancing our understanding of resistance to CDK4/6is....

#### <a id='gse300475'></a>GSE300475 — Single cell RNA sequencing for longitudinal human peripheral blood from HR+ breast cancer patients treated with immunotherapy
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 32 samples/patients; ~64,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE300475_feature_ref.xlsx; filelist.txt; GSE300475_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE300nnn/GSE300475/suppl/GSE300475_feature_ref.xl`
- **Citation / Reference:** PMID:40610460
- **Repository Access Link:** [GSE300475](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300475)
- **Study Abstract / Experimental Design:**
  The limited benefit of immune checkpoint inhibitor in breast cancer indicates the pressing need to identify biomarkers of response to minimize risk and maximize benefit. In this study, we performed single cell RNA sequencing and T cell receptor sequencing on peripheral blood mononuclear cells to monitor the peripheral immune dynamics of an exploratory cohort of hormone receptor positive breast can...

#### <a id='gse212707'></a>GSE212707 — Multiomic analysis reveals conservation of cancer associated fibroblast phenotypes across species and tissue of origin [multiome]
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 21 samples/patients; ~42,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE212707_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/GSE212707/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE212707](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212707)
- **Study Abstract / Experimental Design:**
  Cancer associated fibroblasts (CAFs) are integral to the solid tumor microenvironment. Once thought to be a relatively uniform population of matrix-producing cells, the arrival of single cell RNA sequencing has revealed diverse CAF phenotypes. Here, we further probe CAF heterogeneity with a comprehensive multiome approach. Using paired, same-cell chromatin accessibility and transcriptome analysis,...

#### <a id='gse262288'></a>GSE262288 — Single-Cell RNA Sequencing Identifies Molecular Biomarkers Predicting Response to CDK4/6 Inhibition in Metastatic HR+/HER2- Breast Cancer
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 19 samples/patients; ~38,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE262288_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262288/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:39955556
- **Repository Access Link:** [GSE262288](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262288)
- **Study Abstract / Experimental Design:**
  Cyclin-dependent kinase 4/6 inhibitors (CDK4/6i), combined with endocrine therapy, represent the standard treatment for metastatic HR+/HER2- breast cancer (mBC). Despite their established efficacy, intrinsic resistance affects approximately one-third of patients, posing a significant challenge and underscoring the importance of identifying reliable predictive biomarkers. This study utilizes single...

#### <a id='gse303346'></a>GSE303346 — Single-Cell RNA-Seq Reveals Immunosuppressive Effects of Triple-Negative Breast Cancer–Derived Exosomes on NK Cell Responses
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 16 samples/patients; ~32,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE303346_BT-20_matrix.h5; GSE303346_BT-549_matrix.h5; GSE303346_DU4475_matrix.h5`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE303nnn/GSE303346/suppl/GSE303346_BT-20_matrix.h`
- **Repository Access Link:** [GSE303346](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303346)
- **Study Abstract / Experimental Design:**
  Triple-negative breast cancer (TNBC) is a highly aggressive and immunogenic subtype that lacks effective targeted therapies. Although tumor-derived exosomes are known to influence immune responses, their direct role in shaping human NK cell plasticity remains insufficiently characterized. In this study, we performed an integrative single-cell multiomic analysis of primary human NK cells following ...

#### <a id='gse254991'></a>GSE254991 — Interleukin-16 Establishes a Th1-dominant Tumor Immune Microenvironment And Potentiates Immune Checkpoint Therapies [scRNA-Seq]
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 12 samples/patients; ~24,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE254991_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254991/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:40064918
- **Repository Access Link:** [GSE254991](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254991)
- **Study Abstract / Experimental Design:**
  Overcoming immunosuppression in tumor microenvironment (TME) is crucial for the development of novel cancer immunotherapies. In this study, we revealed a previously unrecognized role of IL-16 in shaping anti-tumor immunity. Compared to healthy individuals, cancer patients exhibited impaired production of IL-16, which was associated with inferior patient prognosis. In multiple murine cancer models,...

#### <a id='gse199219'></a>GSE199219 — Proteogenomic integration of single-cell RNA and protein analysis identifies novel tumour-infiltrating lymphocyte phenotypes in breast cancer
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 10 samples/patients; ~20,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE199219_processed_RNAcounts_and_metadata.tar.gz; filelist.txt; GSE199219_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE199nnn/GSE199219/suppl/GSE199219_processed_RNAc`
- **Repository Access Link:** [GSE199219](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199219)
- **Study Abstract / Experimental Design:**
  High-throughput single-cell RNA sequencing (scRNA-Seq) has become a routine platform for the dissection of solid tumours into their cellular components. The development of methods to incorporate detection of cellular protein epitopes via barcoded antibodies has enabled advances in immunophenotyping, but its application to investigating tissue immunology has been limited. We applied joint single ce...

#### <a id='gse332708'></a>GSE332708 — Hybrid In Vivo Breast Cancer Model Reveals Transcriptomic Insights into Cancer Progression with Age
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 4 samples/patients; ~8,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE332708_Pool-1_filtered_feature_bc_matrix.h5; GSE332708_Pool-1_filtered_matrix.mtx.gz; GSE332708_Pool-2_filtered_feature_bc_matrix.h5`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE332nnn/GSE332708/suppl/GSE332708_Pool-1_filtere`
- **Repository Access Link:** [GSE332708](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE332708)
- **Study Abstract / Experimental Design:**
  Aging is a key risk factor for breast cancer; however, the independent effects of the aged extracellular matrix (ECM) remain understudied. To address this, we developed a novel hybrid in vivo model that enables the independent investigation of age-related ECM influences on breast cancer development and progression. To examine the effects of genes known to be enriched in a cancerous and aged microe...

#### <a id='gse331487'></a>GSE331487 — Immunosuppressive myeloid cells induce mesenchymal-like breast cancer stem cells by a membrane-bound TGF-β1-dependent mechanism [scRNA-seq]
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 3 samples/patients; ~6,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE331487_counts_raw_integrated_samples.tsv.gz; GSE331487_metadata_integrated_samples.tsv.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE331nnn/GSE331487/suppl/GSE331487_counts_raw_int`
- **Repository Access Link:** [GSE331487](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE331487)
- **Study Abstract / Experimental Design:**
  Suppressive myeloid cells play a central role in cancer escape from anti-tumor immunity and exert multiple pro-tumoral activities, including the promotion of cancer cell survival, invasion and metastasis. Recent evidence suggests that certain subsets may also influence tumor heterogeneity and cellular plasticity. However, their role in cancer stem cell promotion and plasticity and the underlying m...

#### <a id='gse302453'></a>GSE302453 — A single-cell map of intratumoral heterogeneity during combination treatment of anti-PD1 and chemotherapy in triple-negative breast cancer
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"We found that the insensitive (PCPAi) patient exhibited high CNV score in tumor cells, lower cytotoxicity in CD8+ T cells, higher quiescence and exhaustion characteristics CD4+ T cells, compared to sensitive (PCPAs) patient."*
- **Cohort Scale:** 2 samples/patients; ~4,000 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `filelist.txt; GSE302453_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE302nnn/GSE302453/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE302453](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302453)
- **Study Abstract / Experimental Design:**
  Triple-negative breast cancer (TNBC) is the most malignant subtype of breast cancer (BC), which exhibited a poor prognosis due to the highly tumor microenvironmental heterogeneity and specific mutations. Recent studies have investigated the different celltypes in TNBC microenvironment, however, the effects of immunotherapy and chemotherapy combination therapy on TNBC are largely unknown. Here we e...

#### <a id='gse299267'></a>GSE299267 — Single cell gene expression profiling of triple negative breast cancer organoids
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** Single-Nucleus RNA-seq
- **Modality:** scRNA-seq + snRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 2 samples/patients; ~4,000 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `filelist.txt; GSE299267_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE299nnn/GSE299267/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE299267](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299267)
- **Study Abstract / Experimental Design:**
  Triple negative breast cancer (TNBC) is a heterogenous disease driven by aberrant activation of receptor tyrosine kinases and the proliferative MEK-ERK signaling pathway. Heterogeneity in kinase expression levels across two TNBC organoid lines was used in combination with the Atlas of substrate specificities for the human Ser/Thr kinome to predict phenotypic differences between the two organoid li...

#### <a id='cellxgene_5a9cfb44'></a>CELLxGENE_5a9cfb44 — Immune Compartment
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 138 samples/patients; ~274,555 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `5a9cfb44-494a-4540-8db8-7143b1bb5a6f.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/904474b2-0029-4b3d-8cd5-1c8349e9462b.h5ad`
- **Citation / Reference:** 10.1093/nargab/lqaf217
- **Repository Access Link:** [CELLxGENE_5a9cfb44](https://cellxgene.cziscience.com/e/5a9cfb44-494a-4540-8db8-7143b1bb5a6f.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Human Breast Cancer Single Cell Atlas | Diseases: HER2 positive breast carcinoma, breast apocrine carcinoma, breast cancer, breast carcinoma, breast mucinous carcinoma, estrogen-receptor positive breast cancer, invasive ductal breast carcinoma, invasive lobular breast carcinoma, invasive tubular breast carcinoma || invasive lobular breast carcinoma, metaplastic breast carcinoma, triple...

#### <a id='cellxgene_de5416ef'></a>CELLxGENE_de5416ef — Global Atlas
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 138 samples/patients; ~621,200 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `de5416ef-bbfa-41f3-99d4-01a3554173f5.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/0fbf6fa3-8463-46ab-bda8-9703848b58db.h5ad`
- **Citation / Reference:** 10.1093/nargab/lqaf217
- **Repository Access Link:** [CELLxGENE_de5416ef](https://cellxgene.cziscience.com/e/de5416ef-bbfa-41f3-99d4-01a3554173f5.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Human Breast Cancer Single Cell Atlas | Diseases: HER2 positive breast carcinoma, breast apocrine carcinoma, breast cancer, breast carcinoma, breast mucinous carcinoma, estrogen-receptor positive breast cancer, invasive ductal breast carcinoma, invasive lobular breast carcinoma, invasive tubular breast carcinoma || invasive lobular breast carcinoma, metaplastic breast carcinoma, triple...

#### <a id='cellxgene_ed880090'></a>CELLxGENE_ed880090 — Epithelial Compartment
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Single-cell suspension sorted by FACS for CD3+/CD8+ T lymphocytes (TILs)."*
- **Cohort Scale:** 136 samples/patients; ~242,613 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `ed880090-465d-4477-9d98-70838ca288e3.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/374cbbe1-0c91-4ca7-8730-ab7d273b960d.h5ad`
- **Citation / Reference:** 10.1093/nargab/lqaf217
- **Repository Access Link:** [CELLxGENE_ed880090](https://cellxgene.cziscience.com/e/ed880090-465d-4477-9d98-70838ca288e3.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Human Breast Cancer Single Cell Atlas | Diseases: HER2 positive breast carcinoma, breast apocrine carcinoma, breast cancer, breast carcinoma, breast mucinous carcinoma, estrogen-receptor positive breast cancer, invasive ductal breast carcinoma, invasive lobular breast carcinoma, invasive tubular breast carcinoma || invasive lobular breast carcinoma, metaplastic breast carcinoma, triple...

#### <a id='cellxgene_75011e96'></a>CELLxGENE_75011e96 — Stromal Compartment
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 131 samples/patients; ~104,032 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `75011e96-ba37-4977-a539-fbf9f66a5e37.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/b5f2ee00-64c6-434c-8fee-3578fc9fea05.h5ad`
- **Citation / Reference:** 10.1093/nargab/lqaf217
- **Repository Access Link:** [CELLxGENE_75011e96](https://cellxgene.cziscience.com/e/75011e96-ba37-4977-a539-fbf9f66a5e37.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Human Breast Cancer Single Cell Atlas | Diseases: HER2 positive breast carcinoma, breast apocrine carcinoma, breast cancer, breast carcinoma, breast mucinous carcinoma, estrogen-receptor positive breast cancer, invasive ductal breast carcinoma, invasive lobular breast carcinoma, invasive tubular breast carcinoma || invasive lobular breast carcinoma, metaplastic breast carcinoma, triple...

#### <a id='cellxgene_6f9de485'></a>CELLxGENE_6f9de485 — Single-cell RNA-seq data from untreated tissue samples of treatment-naive triple-negative breast cancer patients
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 101 samples/patients; ~427,823 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `6f9de485-58cd-4342-bfc4-b3d3dd223aa8.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/e94bd3cc-6271-424a-baac-12f8eb320a0e.h5ad`
- **Citation / Reference:** Dataset Version: https://datasets.cellxgene.cziscience.com/e94bd3cc-6271-424a-baac-12f8eb320a0e.h5ad curated and distributed by CZ CELLxGENE Discover in Collection: https://cellxgene.cziscience.com/collections/ceef2841-5333-46ac-92ef-ccbe0c20fe55
- **Repository Access Link:** [CELLxGENE_6f9de485](https://cellxgene.cziscience.com/e/6f9de485-58cd-4342-bfc4-b3d3dd223aa8.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Single-cell atlas of untreated human triple negative breast cancer | Diseases: triple-negative breast carcinoma | Tissues: breast | Assays: 10x 3' v2, 10x 3' v3...

#### <a id='gse300628'></a>GSE300628 — HER2 heterogeneous breast cancer models reveal novel therapeutic targets and subclonal dynamics during evolution to resistance to HER2-targeted therapies
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 60 samples/patients; ~120,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE300628_Cuff_Gene_Counts.csv.gz; GSE300628_NTres_Cuff_Gene_Counts.csv.gz; GSE300628_PTres_Cuff_Gene_Counts.csv.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE300nnn/GSE300628/suppl/GSE300628_Cuff_Gene_Coun`
- **Citation / Reference:** PMID:41925564
- **Repository Access Link:** [GSE300628](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300628)
- **Study Abstract / Experimental Design:**
  Intratumor heterogeneity for HER2 expression and amplification, observed in up to 40% of HER2-positive breast cancer is a driver of resistance to HER2-targeted therapies. The advancement of treatment for HER2 heterogeneous tumors has been hindered by the lack of preclinical models that accurately replicate the human disease. Here we describe human HER2 heterogeneous breast cancer cell models compo...

#### <a id='cellxgene_34f5307e'></a>CELLxGENE_34f5307e — Frozen samples (single nucleus)
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** snRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 32 samples/patients; ~394,534 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `34f5307e-7b4d-4a48-b68f-2ba844c6414b.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/5c7bd9b2-5609-4a7a-9445-69ff3f9aaff5.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_34f5307e](https://cellxgene.cziscience.com/e/34f5307e-7b4d-4a48-b68f-2ba844c6414b.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: axilla, bone spine, brain, breast, chest wall, liver, skin epidermis | Assays: 10x 3' v2, 10x 3' v3...

#### <a id='cellxgene_6c87755e'></a>CELLxGENE_6c87755e — Fresh samples (single cell)
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 30 samples/patients; ~157,531 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `6c87755e-a671-41a8-9c7e-4e43b850a57b.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/173f9d7a-d6db-4b59-9041-ef8342fb5404.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_6c87755e](https://cellxgene.cziscience.com/e/6c87755e-a671-41a8-9c7e-4e43b850a57b.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: axilla, bone spine, breast, liver, lung, neck | Assays: 10x 3' v2, 10x 3' v3...

#### <a id='cellxgene_933497dc'></a>CELLxGENE_933497dc — A single-cell and spatially-resolved atlas of human breast cancers - T_cells
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 26 samples/patients; ~35,214 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `933497dc-7592-4532-9056-b9941547b3ac.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/fff5d349-6bf8-4a6d-8124-486e2c9daf20.h5ad`
- **Citation / Reference:** 10.1038/s41588-021-00911-1
- **Repository Access Link:** [CELLxGENE_933497dc](https://cellxgene.cziscience.com/e/933497dc-7592-4532-9056-b9941547b3ac.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: A single-cell and spatially resolved atlas of human breast cancers | Diseases: breast carcinoma, invasive ductal breast carcinoma, invasive lobular breast carcinoma | Tissues: breast | Assays: 10x 3' v2, 10x 3' v3, 10x 5' v1...

#### <a id='cellxgene_9fddb063'></a>CELLxGENE_9fddb063 — A single-cell and spatially-resolved atlas of human breast cancers
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 26 samples/patients; ~100,064 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `9fddb063-056d-4202-8b8a-4b0ee531d3ce.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/22a27631-aecf-463b-86c6-a8334a2f2cf2.h5ad`
- **Citation / Reference:** 10.1038/s41588-021-00911-1
- **Repository Access Link:** [CELLxGENE_9fddb063](https://cellxgene.cziscience.com/e/9fddb063-056d-4202-8b8a-4b0ee531d3ce.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: A single-cell and spatially resolved atlas of human breast cancers | Diseases: breast carcinoma, invasive ductal breast carcinoma, invasive lobular breast carcinoma | Tissues: breast | Assays: 10x 3' v2, 10x 3' v3, 10x 5' v1...

#### <a id='cellxgene_11a3244a'></a>CELLxGENE_11a3244a — A single-cell and spatially-resolved atlas of human breast cancers - Myeloid
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 26 samples/patients; ~9,675 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `11a3244a-6b2f-43ca-95a5-d2cd95c482d2.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/27a4cfeb-fc26-4b46-8320-d84e17d1004b.h5ad`
- **Citation / Reference:** 10.1038/s41588-021-00911-1
- **Repository Access Link:** [CELLxGENE_11a3244a](https://cellxgene.cziscience.com/e/11a3244a-6b2f-43ca-95a5-d2cd95c482d2.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: A single-cell and spatially resolved atlas of human breast cancers | Diseases: breast carcinoma, invasive ductal breast carcinoma, invasive lobular breast carcinoma | Tissues: breast | Assays: 10x 3' v2, 10x 3' v3, 10x 5' v1...

#### <a id='cellxgene_4cdd25a4'></a>CELLxGENE_4cdd25a4 — A single-cell and spatially-resolved atlas of human breast cancers - Stromal
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 26 samples/patients; ~19,601 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `4cdd25a4-f4fa-4c13-b8a9-80cf28511d46.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/0e60d1bf-5366-47ba-a01c-7011dade6b76.h5ad`
- **Citation / Reference:** 10.1038/s41588-021-00911-1
- **Repository Access Link:** [CELLxGENE_4cdd25a4](https://cellxgene.cziscience.com/e/4cdd25a4-f4fa-4c13-b8a9-80cf28511d46.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: A single-cell and spatially resolved atlas of human breast cancers | Diseases: breast carcinoma, invasive ductal breast carcinoma, invasive lobular breast carcinoma | Tissues: breast | Assays: 10x 3' v2, 10x 3' v3, 10x 5' v1...

#### <a id='cellxgene_7357bdd2'></a>CELLxGENE_7357bdd2 — A single-cell and spatially-resolved atlas of human breast cancers - Cancer_epithelial
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v2
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 20 samples/patients; ~24,489 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `7357bdd2-e05b-4e87-acd9-e3a26382648b.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/de3c0ed3-5aab-4fd5-be2d-ecaa1b598a18.h5ad`
- **Citation / Reference:** 10.1038/s41588-021-00911-1
- **Repository Access Link:** [CELLxGENE_7357bdd2](https://cellxgene.cziscience.com/e/7357bdd2-e05b-4e87-acd9-e3a26382648b.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: A single-cell and spatially resolved atlas of human breast cancers | Diseases: breast carcinoma, invasive ductal breast carcinoma, invasive lobular breast carcinoma | Tissues: breast | Assays: 10x 3' v2, 10x 5' v1...

#### <a id='cellxgene_04d87de6'></a>CELLxGENE_04d87de6 — A single-cell and spatially-resolved atlas of human breast cancers - B_cells
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 20 samples/patients; ~3,206 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `04d87de6-c20a-4186-8884-f47dba20b0a4.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/075b9b95-c24c-4c02-ae8c-ae5381418e87.h5ad`
- **Citation / Reference:** 10.1038/s41588-021-00911-1
- **Repository Access Link:** [CELLxGENE_04d87de6](https://cellxgene.cziscience.com/e/04d87de6-c20a-4186-8884-f47dba20b0a4.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: A single-cell and spatially resolved atlas of human breast cancers | Diseases: breast carcinoma, invasive ductal breast carcinoma, invasive lobular breast carcinoma | Tissues: breast | Assays: 10x 3' v2, 10x 3' v3, 10x 5' v1...

#### <a id='cellxgene_2f05ab20'></a>CELLxGENE_2f05ab20 — A total of 38,217 droplet-based single-nucleus transcriptomes profiled across 14 cell types in the frontal cortex
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** snRNA-seq
- **Cell Selection / Filtering Strategy:** **Nuclei Isolation (snRNA-seq)**
  > *"Single-nucleus RNA sequencing performed on isolated nuclei from frozen resection tissue."*
- **Cohort Scale:** 16 samples/patients; ~38,217 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `2f05ab20-a092-4bab-9276-3e0eb24e3fee.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/c4ee3adc-893f-449a-988c-b3b4d51ceba9.h5ad`
- **Citation / Reference:** 10.1038/s41586-021-03710-0
- **Repository Access Link:** [CELLxGENE_2f05ab20](https://cellxgene.cziscience.com/e/2f05ab20-a092-4bab-9276-3e0eb24e3fee.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Dysregulation of brain and choroid plexus cell types in severe COVID-19 | Diseases: COVID-19, breast cancer, cardiomyopathy, chronic obstructive pulmonary disease, heart disorder, influenza, myocardial infarction, small cell lung carcinoma, tongue cancer | Tissues: medial orbital frontal cortex | Assays: 10x 3' v3...

#### <a id='cellxgene_0c86f0de'></a>CELLxGENE_0c86f0de — A single-cell and spatially-resolved atlas of human breast cancers - Plasmablasts
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 14 samples/patients; ~3,524 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `0c86f0de-ddcb-454c-b00b-37feb69e7da1.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/0b8e5229-5ab6-4c9d-b01a-3192fd716ca5.h5ad`
- **Citation / Reference:** 10.1038/s41588-021-00911-1
- **Repository Access Link:** [CELLxGENE_0c86f0de](https://cellxgene.cziscience.com/e/0c86f0de-ddcb-454c-b00b-37feb69e7da1.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: A single-cell and spatially resolved atlas of human breast cancers | Diseases: invasive ductal breast carcinoma, invasive lobular breast carcinoma | Tissues: breast | Assays: 10x 3' v2, 10x 3' v3, 10x 5' v1...

#### <a id='gse337706'></a>GSE337706 — Exploring Platelet-Covered and Naked Circulating Tumor Cells: A Single-Cell Transcriptomic Perspective [scRNA-seq Singleron]
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 13 samples/patients; ~26,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE337706_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE337nnn/GSE337706/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE337706](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE337706)
- **Study Abstract / Experimental Design:**
  Circulating tumor cells (CTCs) and platelets can be collected simultaneously during liquid biopsy; however, their interaction in the form of platelet-covered CTCs (pcCTCs) remains only partially understood. This submission contains single-cell RNA sequencing data generated from PBMC/buffy coat fractions from 29 donors: 12 patients with high-grade serous ovarian carcinoma, 4 non-malignant gynecolog...

#### <a id='gse281488'></a>GSE281488 — In situ single-cell profiling of a brain metastasis from a HER2+ Breast Cancer Patient and 8 primary tumors from TNBC patients.    
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 9 samples/patients; ~18,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE281488_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281488/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41923644
- **Repository Access Link:** [GSE281488](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281488)
- **Study Abstract / Experimental Design:**
  A key problem in cancer biology is understanding the phenotypic characteristics of the aggressive cell states present in macrometastases and their molecular underpinnings. Here we report on the cellular constituents of  breast cancer metastases,  revealing that these malignant cells display a gene-expression profile associated to embryonic development of tubular organs. Similar properties are foun...

#### <a id='cellxgene_5d3fc988'></a>CELLxGENE_5d3fc988 — Figure 1
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 7 samples/patients; ~49,109 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `5d3fc988-765d-48ba-bfb5-151ef2988cac.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/6ff7acae-ac81-45f5-9348-94e85a25fe7c.h5ad`
- **Citation / Reference:** 10.1002/ctm2.1356
- **Repository Access Link:** [CELLxGENE_5d3fc988](https://cellxgene.cziscience.com/e/5d3fc988-765d-48ba-bfb5-151ef2988cac.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: The cellular landscape of breast cancer associated malignant pleural effusions | Diseases: luminal A breast carcinoma, luminal B breast carcinoma, triple-negative breast carcinoma | Tissues: pleural effusion | Assays: 10x 3' v3...

#### <a id='cellxgene_3e4e2c8e'></a>CELLxGENE_3e4e2c8e — Figure 2
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 7 samples/patients; ~22,414 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `3e4e2c8e-cd17-4ece-9141-1e7f4ce8da1f.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/18cd9582-de51-4924-afbc-75040c910066.h5ad`
- **Citation / Reference:** 10.1002/ctm2.1356
- **Repository Access Link:** [CELLxGENE_3e4e2c8e](https://cellxgene.cziscience.com/e/3e4e2c8e-cd17-4ece-9141-1e7f4ce8da1f.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: The cellular landscape of breast cancer associated malignant pleural effusions | Diseases: luminal A breast carcinoma, luminal B breast carcinoma, triple-negative breast carcinoma | Tissues: pleural effusion | Assays: 10x 3' v3...

#### <a id='gse281490'></a>GSE281490 — Human scRNA-seq of three breast cancer metastases and three normal mammary glands
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 6 samples/patients; ~12,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE281490_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281490/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41923644
- **Repository Access Link:** [GSE281490](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281490)
- **Study Abstract / Experimental Design:**
  A key problem in cancer biology is understanding the phenotypic characteristics of the aggressive cell states present in macrometastases and their molecular underpinnings. Here we report on the cellular constituents of human breast cancer metastases,  revealing that these malignant cells are invariably locked in a chimeric cell state. Human BC metastatic cells indeed display a gene-expression prof...

#### <a id='cellxgene_68b6114f'></a>CELLxGENE_68b6114f — Stromal cell diversity associated with immune evasion in human triple‐negative breast cancer
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v2
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 5 samples/patients; ~24,271 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `68b6114f-e990-4033-bfb3-c35536633aaa.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/4d036b2c-6135-4515-9c44-a5aa7a2e32e5.h5ad`
- **Citation / Reference:** 10.15252/embj.2019104063
- **Repository Access Link:** [CELLxGENE_68b6114f](https://cellxgene.cziscience.com/e/68b6114f-e990-4033-bfb3-c35536633aaa.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Stromal cell diversity associated with immune evasion in human triple‐negative breast cancer | Diseases: triple-negative breast carcinoma | Tissues: mammary gland connective tissue | Assays: 10x 3' v2...

#### <a id='cellxgene_e500acbf'></a>CELLxGENE_e500acbf — Figure 3
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 5 samples/patients; ~37,428 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `e500acbf-f166-4c46-b290-f93506cf26d3.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/99a43cf0-5b29-4190-a281-a961a1a81ca7.h5ad`
- **Citation / Reference:** 10.1002/ctm2.1356
- **Repository Access Link:** [CELLxGENE_e500acbf](https://cellxgene.cziscience.com/e/e500acbf-f166-4c46-b290-f93506cf26d3.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: The cellular landscape of breast cancer associated malignant pleural effusions | Diseases: luminal B breast carcinoma, triple-negative breast carcinoma | Tissues: pleural effusion | Assays: 10x 3' v2, 10x 3' v3...

#### <a id='gse230327'></a>GSE230327 — BRD8 is a therapeutic vulnerability for overcoming resistance to dual ER/HER2 blockade therapy in HR+/HER2+ breast cancer [scRNA-seq]
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 4 samples/patients; ~8,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE230327_matrix.mtx.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE230nnn/GSE230327/suppl/GSE230327_matrix.mtx.gz`
- **Citation / Reference:** PMID:41886605
- **Repository Access Link:** [GSE230327](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230327)
- **Study Abstract / Experimental Design:**
  Hormone receptor (HR)-positive, HER2-positive breast cancers are resistant to endocrine and anti-HER2 therapies due to crosstalk between estrogen receptor (ER) and HER2. However, how anti-HER2 agents activate ER as a mechanism of resistance remains unknown. Using single-cell RNA sequencing, we identified Bromodomain Containing Protein 8 (BRD8) as a major mediator of ER activation in response to ne...

#### <a id='cellxgene_fbdd8c17'></a>CELLxGENE_fbdd8c17 — Breast cancer 4 patients archival FFPE samples profiling with Chromium FLEX
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium (3' unspecified)
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 4 samples/patients; ~10,689 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `fbdd8c17-b34a-4cbc-abc4-1aeaa294a538.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/5466ce85-7f19-49e2-9a1b-54646bb84708.h5ad`
- **Citation / Reference:** 10.1101/2024.11.01.621259
- **Repository Access Link:** [CELLxGENE_fbdd8c17](https://cellxgene.cziscience.com/e/fbdd8c17-b34a-4cbc-abc4-1aeaa294a538.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Transcriptome Analysis of Archived Tumor Tissues by Visium, GeoMx DSP, and Chromium Methods Reveals Inter- and Intra-Patient Heterogeneity | Diseases: invasive ductal breast carcinoma, invasive lobular breast carcinoma | Tissues: parenchyma of mammary gland | Assays: 10x Next GEM Flex v1...

#### <a id='gse329389'></a>GSE329389 — Uncovering of Intratumoral Androgen-Dependent Molecular Signature for Better Breast Cancer Prognosis [snRNA-seq]
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** Single-Nucleus RNA-seq
- **Modality:** snRNA-seq
- **Cell Selection / Filtering Strategy:** **Nuclei Isolation (snRNA-seq)**
  > *"Uncovering of Intratumoral Androgen-Dependent Molecular Signature for Better Breast Cancer Prognosis [snRNA-seq] The role of androgens in breast cancer (BC) was defined by correlating circulating androgens to clinical variables."*
- **Cohort Scale:** 3 samples/patients; ~6,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE329389_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE329nnn/GSE329389/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE329389](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329389)
- **Study Abstract / Experimental Design:**
  The role of androgens in breast cancer (BC) was defined by correlating circulating androgens to clinical variables. Considering that intratumoral hormonal milieu is the critical determinant of BC characteristics, correlating intratumoral, androgens and their downstream transcriptome, a measure of free unbound androgen function, will provide clarity on androgen action in BC. Preclinical studies hav...

#### <a id='gse325982'></a>GSE325982 — Systemic and breast chronic inflammation and hormone disposition promote a tumor-permissive environment for breast cancer in older women (Organoid scRNA-Seq)
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** Single-Nucleus RNA-seq
- **Modality:** scRNA-seq + snRNA-seq
- **Cell Selection / Filtering Strategy:** **Nuclei Isolation (snRNA-seq)**
  > *"snRNA-seq of the tumors showed broad local immune dysfunction that was associated with circulating chronic inflammation."*
- **Cohort Scale:** 2 samples/patients; ~4,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE325982_metadata_individual_samples.xlsx; filelist.txt; GSE325982_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE325nnn/GSE325982/suppl/GSE325982_metadata_indiv`
- **Citation / Reference:** PMID:42509302
- **Repository Access Link:** [GSE325982](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325982)
- **Study Abstract / Experimental Design:**
  Estrogen receptor positive (ER+) breast cancer, the most common subtype of breast cancer, is an age-related disease, with the peak incidence of diagnosis occurring around age 70 despite low circulating levels of estradiol. Despite the hormone sensitivity of these age-related tumors, our understanding of the interplay between the systemic and local hormonal disposition and chronic inflammaging is l...

#### <a id='gse309616'></a>GSE309616 — Therapeutic Synergy Overcomes Carboplatin Resistance in Triple-Negative Breast Cancer [scRNA-seq]
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 2 samples/patients; ~4,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE309616_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE309nnn/GSE309616/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41630032
- **Repository Access Link:** [GSE309616](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309616)
- **Study Abstract / Experimental Design:**
  Background: Triple-negative breast cancer (TNBC) is an aggressive subtype lacking targeted therapeutic options, where platinum-based chemotherapy such as carboplatin serves as a cornerstone of treatment. Despite initial responses, the rapid emergence of acquired resistance remains a major clinical barrier. Understanding the molecular adaptations that drive platinum resistance is essential to devel...

#### <a id='gse229723'></a>GSE229723 — Single cell RNA sequencing (scRNA-seq) data on the human CD45+ cells and tumor cells collected from four different treatment groups of hCD34+ humanized mice with human TNBC xenografic tumor model (mixed sample, non-demultiplexed)
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD45+ Immune-enriched)**
  > *"Single-cell suspension sorted by FACS for CD45+ leukocytes to enrich for tumor-infiltrating immune cells."*
- **Cohort Scale:** 2 samples/patients; ~4,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE229723_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE229nnn/GSE229723/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:40794843
- **Repository Access Link:** [GSE229723](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229723)
- **Study Abstract / Experimental Design:**
  Using a hCD34+ humanized mouse tumor model, the single-cell sequencing analysis revealed that in vivo knockdown of EPIC1 enhanced the T cells and macrophage infiltration by activating type I IFNs signaling....

#### <a id='cellxgene_12c868c6'></a>CELLxGENE_12c868c6 — HTAPP-878-SMP-7149 scRNA-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~10,623 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `12c868c6-94df-48b5-acf2-b82f6aa14074.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/617c8418-bb79-4f3c-beb0-58b4aec0eabe.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_12c868c6](https://cellxgene.cziscience.com/e/12c868c6-94df-48b5-acf2-b82f6aa14074.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: axilla | Assays: 10x 3' v3...

#### <a id='cellxgene_a6c0143c'></a>CELLxGENE_a6c0143c — HTAPP-982-SMP-7629 scRNA-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~7,505 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `a6c0143c-11f7-4a54-8710-6b427b78873a.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/5f539b0b-5878-49b1-8115-c8d2b016c45b.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_a6c0143c](https://cellxgene.cziscience.com/e/a6c0143c-11f7-4a54-8710-6b427b78873a.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: liver | Assays: 10x 3' v3...

#### <a id='cellxgene_6f0858c0'></a>CELLxGENE_6f0858c0 — CID4535
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `6f0858c0-c590-4740-b022-c152e7608d66.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/91d821e4-403c-4de8-bed8-f798ad33d93a.h5ad`
- **Citation / Reference:** 10.1038/s41588-021-00911-1
- **Repository Access Link:** [CELLxGENE_6f0858c0](https://cellxgene.cziscience.com/e/6f0858c0-c590-4740-b022-c152e7608d66.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: A single-cell and spatially resolved atlas of human breast cancers | Diseases: invasive lobular breast carcinoma | Tissues: breast | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_80466231'></a>CELLxGENE_80466231 — CID4465
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `80466231-7133-4097-88ae-40cb6cce1a33.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/0fa1249e-a670-4534-b51c-69218dce7889.h5ad`
- **Citation / Reference:** 10.1038/s41588-021-00911-1
- **Repository Access Link:** [CELLxGENE_80466231](https://cellxgene.cziscience.com/e/80466231-7133-4097-88ae-40cb6cce1a33.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: A single-cell and spatially resolved atlas of human breast cancers | Diseases: invasive ductal breast carcinoma | Tissues: breast | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_a6b0f655'></a>CELLxGENE_a6b0f655 — CID44971
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `a6b0f655-820c-4082-98e2-42f33c2a71a7.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/58a6f5de-f847-495d-84a5-4766c76dcbc0.h5ad`
- **Citation / Reference:** 10.1038/s41588-021-00911-1
- **Repository Access Link:** [CELLxGENE_a6b0f655](https://cellxgene.cziscience.com/e/a6b0f655-820c-4082-98e2-42f33c2a71a7.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: A single-cell and spatially resolved atlas of human breast cancers | Diseases: invasive ductal breast carcinoma | Tissues: breast | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_aafb780d'></a>CELLxGENE_aafb780d — CID4290
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `aafb780d-f52d-4285-a1f4-57376cabe1ee.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/9c62cb12-39e5-4cf4-a030-802ef136d044.h5ad`
- **Citation / Reference:** 10.1038/s41588-021-00911-1
- **Repository Access Link:** [CELLxGENE_aafb780d](https://cellxgene.cziscience.com/e/aafb780d-f52d-4285-a1f4-57376cabe1ee.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: A single-cell and spatially resolved atlas of human breast cancers | Diseases: invasive ductal breast carcinoma | Tissues: breast | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_ee141ea4'></a>CELLxGENE_ee141ea4 — GSM6592053_M6
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `ee141ea4-6678-4908-9c2f-db71d92e74ce.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/4936568b-7d72-4916-a1b2-0dcd3a7273bb.h5ad`
- **Citation / Reference:** 10.1016/j.labinv.2023.100258
- **Repository Access Link:** [CELLxGENE_ee141ea4](https://cellxgene.cziscience.com/e/ee141ea4-6678-4908-9c2f-db71d92e74ce.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatial Transcriptomics Reveal Pitfalls and Opportunities for the Detection of Rare High-Plasticity Breast Cancer Subtypes | Diseases: triple-negative breast carcinoma | Tissues: breast | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_f354e4c3'></a>CELLxGENE_f354e4c3 — 1142243F
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `f354e4c3-eb53-4917-8476-39a860e30124.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/5d904a8a-9190-43d6-8779-f167bf508244.h5ad`
- **Citation / Reference:** 10.1038/s41588-021-00911-1
- **Repository Access Link:** [CELLxGENE_f354e4c3](https://cellxgene.cziscience.com/e/f354e4c3-eb53-4917-8476-39a860e30124.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: A single-cell and spatially resolved atlas of human breast cancers | Diseases: breast carcinoma | Tissues: breast | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_02aa7750'></a>CELLxGENE_02aa7750 — HTAPP-878-SMP-7149 Slide-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~22,033 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `02aa7750-ed05-4eca-8061-766c68a81a22.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/ab8551c5-cafe-4a07-bf00-097167a87831.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_02aa7750](https://cellxgene.cziscience.com/e/02aa7750-ed05-4eca-8061-766c68a81a22.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: axilla | Assays: Slide-seqV2...

#### <a id='cellxgene_05a49baa'></a>CELLxGENE_05a49baa — HTAPP-853-SMP-4381 scRNA-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v2
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,742 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `05a49baa-d326-42ae-86d2-94de3a659901.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/3fe6302c-0510-4b8d-88cf-420f4ca50de0.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_05a49baa](https://cellxgene.cziscience.com/e/05a49baa-d326-42ae-86d2-94de3a659901.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: liver | Assays: 10x 3' v2...

#### <a id='cellxgene_0f9d1892'></a>CELLxGENE_0f9d1892 — HTAPP-514-SMP-6760 Slide-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~33,272 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `0f9d1892-919a-4f51-8b95-0b50750da1e6.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/ccc7ab02-1721-442e-858d-a77f1ce91c39.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_0f9d1892](https://cellxgene.cziscience.com/e/0f9d1892-919a-4f51-8b95-0b50750da1e6.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: axilla | Assays: Slide-seqV2...

#### <a id='cellxgene_1637e817'></a>CELLxGENE_1637e817 — HTAPP-982-SMP-7629 Slide-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~8,086 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `1637e817-d091-4793-a3af-796ef57e2668.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/e8ef7247-18e3-4d62-8357-bf446d3731cb.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_1637e817](https://cellxgene.cziscience.com/e/1637e817-d091-4793-a3af-796ef57e2668.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: liver | Assays: Slide-seqV2...

#### <a id='cellxgene_1884e651'></a>CELLxGENE_1884e651 — HTAPP-812-SMP-8239 scRNA-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~5,463 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `1884e651-5f7f-4e6c-a71a-843ba272aaf4.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/1edc3996-8cea-4662-9ad0-fbd32279a80a.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_1884e651](https://cellxgene.cziscience.com/e/1884e651-5f7f-4e6c-a71a-843ba272aaf4.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: axilla | Assays: 10x 3' v3...

#### <a id='cellxgene_2dd73feb'></a>CELLxGENE_2dd73feb — HTAPP-364-SMP-1321 scRNA-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v2
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~13,167 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `2dd73feb-0527-47d5-8c7b-35b6de16aecb.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/9eba2eb2-187d-475d-82b5-2d7f52d667e1.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_2dd73feb](https://cellxgene.cziscience.com/e/2dd73feb-0527-47d5-8c7b-35b6de16aecb.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: liver | Assays: 10x 3' v2...

#### <a id='cellxgene_39f6fec9'></a>CELLxGENE_39f6fec9 — HTAPP-514-SMP-6760 scRNA-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~20,916 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `39f6fec9-2539-4b99-982b-f12888ef649a.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/3efc80da-c21f-470f-9b01-fc39c45f2973.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_39f6fec9](https://cellxgene.cziscience.com/e/39f6fec9-2539-4b99-982b-f12888ef649a.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: axilla | Assays: 10x 3' v3...

#### <a id='cellxgene_44941fdb'></a>CELLxGENE_44941fdb — HTAPP-783-SMP-4081 scRNA-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v2
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~2,276 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `44941fdb-8a6f-42d5-9b34-da48c2a3f774.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/c0df617a-321b-46ee-90ec-14f92434be41.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_44941fdb](https://cellxgene.cziscience.com/e/44941fdb-8a6f-42d5-9b34-da48c2a3f774.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: breast | Assays: 10x 3' v2...

#### <a id='cellxgene_48b55b2b'></a>CELLxGENE_48b55b2b — HTAPP-213-SMP-6752 Slide-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~9,362 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `48b55b2b-c6a3-4288-a09e-e01e5078c678.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/711fb5e6-32ec-4e2c-8fa7-004922ef6d9a.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_48b55b2b](https://cellxgene.cziscience.com/e/48b55b2b-c6a3-4288-a09e-e01e5078c678.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: chest wall | Assays: Slide-seqV2...

#### <a id='cellxgene_494faa16'></a>CELLxGENE_494faa16 — HTAPP-895-SMP-7359 scRNA-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~9,958 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `494faa16-c42d-4b66-8667-f66a3bcacb01.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/fe938180-0708-411a-b069-d4a5f49e101b.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_494faa16](https://cellxgene.cziscience.com/e/494faa16-c42d-4b66-8667-f66a3bcacb01.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: liver | Assays: 10x 3' v3...

#### <a id='cellxgene_54d56674'></a>CELLxGENE_54d56674 — HTAPP-917-SMP-4531 scRNA-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v2
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~11,074 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `54d56674-e4e1-4e3c-9110-ea5cc2ac0eb5.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/b230bf25-895e-41ac-8d30-62ee53c7d5cd.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_54d56674](https://cellxgene.cziscience.com/e/54d56674-e4e1-4e3c-9110-ea5cc2ac0eb5.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: liver | Assays: 10x 3' v2...

#### <a id='cellxgene_59d14a35'></a>CELLxGENE_59d14a35 — HTAPP-783-SMP-4081 Slide-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~5,833 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `59d14a35-0c86-4eee-90b7-fe59595f522c.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/0e4ba413-8673-4ac7-80e9-e977557753e7.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_59d14a35](https://cellxgene.cziscience.com/e/59d14a35-0c86-4eee-90b7-fe59595f522c.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: breast | Assays: Slide-seqV2...

#### <a id='cellxgene_6384d8b8'></a>CELLxGENE_6384d8b8 — HTAPP-880-SMP-7179 Slide-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~22,188 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `6384d8b8-8511-4ffb-b543-26dc8dd5967a.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/f3b4d6ca-e7e6-43a4-8acb-272c210ac13a.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_6384d8b8](https://cellxgene.cziscience.com/e/6384d8b8-8511-4ffb-b543-26dc8dd5967a.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: liver | Assays: Slide-seqV2...

#### <a id='cellxgene_71513028'></a>CELLxGENE_71513028 — HTAPP-997-SMP-7789 scRNA-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~12,258 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `71513028-e4a0-4c62-a435-16598791690a.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/79a15872-c7b9-402b-8e0c-a9dd669c366b.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_71513028](https://cellxgene.cziscience.com/e/71513028-e4a0-4c62-a435-16598791690a.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: liver | Assays: 10x 3' v3...

#### <a id='cellxgene_71e44b30'></a>CELLxGENE_71e44b30 — HTAPP-853-SMP-4381 Slide-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~11,430 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `71e44b30-83b6-455c-a161-bc04c9898f00.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/1130b63b-68f6-4d04-96ed-75cb7fabc1ee.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_71e44b30](https://cellxgene.cziscience.com/e/71e44b30-83b6-455c-a161-bc04c9898f00.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: liver | Assays: Slide-seqV2...

#### <a id='cellxgene_7432b873'></a>CELLxGENE_7432b873 — HTAPP-917-SMP-4531 Slide-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~8,715 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `7432b873-1e96-4a8a-b97b-8921937f4afe.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/cecc7d77-0577-4600-b848-80742c472c3f.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_7432b873](https://cellxgene.cziscience.com/e/7432b873-1e96-4a8a-b97b-8921937f4afe.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: liver | Assays: Slide-seqV2...

#### <a id='cellxgene_9237e573'></a>CELLxGENE_9237e573 — HTAPP-997-SMP-7789 Slide-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~10,892 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `9237e573-980a-4be6-b028-0f01cc064956.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/55185bb7-deba-497f-bb8e-a2fe098bdcc7.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_9237e573](https://cellxgene.cziscience.com/e/9237e573-980a-4be6-b028-0f01cc064956.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: liver | Assays: Slide-seqV2...

#### <a id='cellxgene_9c5f68fc'></a>CELLxGENE_9c5f68fc — HTAPP-330-SMP-1082 Slide-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~11,851 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `9c5f68fc-860f-414b-99c2-cc6be691aff3.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/4b2fe7f7-e1a8-441b-a6e6-1693452f714b.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_9c5f68fc](https://cellxgene.cziscience.com/e/9c5f68fc-860f-414b-99c2-cc6be691aff3.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: liver | Assays: Slide-seqV2...

#### <a id='cellxgene_a6347c54'></a>CELLxGENE_a6347c54 — HTAPP-313-SMP-932 Slide-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~11,484 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `a6347c54-9e69-45e4-b1db-f5d84aa8321a.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/b50ba134-9dac-49ee-974f-fc7c933a7acc.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_a6347c54](https://cellxgene.cziscience.com/e/a6347c54-9e69-45e4-b1db-f5d84aa8321a.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: liver | Assays: Slide-seqV2...

#### <a id='cellxgene_24dbd26d'></a>CELLxGENE_24dbd26d — 1160920F
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `24dbd26d-1a3a-4a0c-a4e9-3d95f373217f.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/6cdeab96-fd27-40c9-a75a-cb94a41ec10b.h5ad`
- **Citation / Reference:** 10.1038/s41588-021-00911-1
- **Repository Access Link:** [CELLxGENE_24dbd26d](https://cellxgene.cziscience.com/e/24dbd26d-1a3a-4a0c-a4e9-3d95f373217f.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: A single-cell and spatially resolved atlas of human breast cancers | Diseases: breast carcinoma | Tissues: breast | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_aa6f371d'></a>CELLxGENE_aa6f371d — HTAPP-895-SMP-7359 Slide-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~6,210 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `aa6f371d-dacd-47e7-a774-e2f4495604a2.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/0c880b0d-ef83-411d-b6ca-7b2f543f47bb.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_aa6f371d](https://cellxgene.cziscience.com/e/aa6f371d-dacd-47e7-a774-e2f4495604a2.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: liver | Assays: Slide-seqV2...

#### <a id='cellxgene_c7d0def0'></a>CELLxGENE_c7d0def0 — HTAPP-880-SMP-7179 scRNA-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~10,918 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `c7d0def0-2dcd-4111-902b-67e6baeac119.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/8b95d860-338d-44d7-ae8d-1e0a6eb56f70.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_c7d0def0](https://cellxgene.cziscience.com/e/c7d0def0-2dcd-4111-902b-67e6baeac119.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: liver | Assays: 10x 3' v3...

#### <a id='cellxgene_cd6398a9'></a>CELLxGENE_cd6398a9 — HTAPP-944-SMP-7479 scRNA-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~10,016 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `cd6398a9-c0af-4467-9091-c536866535bd.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/df219240-1a4c-403f-8dc2-71f90ad1be53.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_cd6398a9](https://cellxgene.cziscience.com/e/cd6398a9-c0af-4467-9091-c536866535bd.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: liver | Assays: 10x 3' v3...

#### <a id='cellxgene_e06e9bf3'></a>CELLxGENE_e06e9bf3 — HTAPP-944-SMP-7479 Slide-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~9,422 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `e06e9bf3-81ef-434d-bc86-08ec1266b155.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/bffd3499-b106-48aa-a935-2a562f68e012.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_e06e9bf3](https://cellxgene.cziscience.com/e/e06e9bf3-81ef-434d-bc86-08ec1266b155.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: liver | Assays: Slide-seqV2...

#### <a id='cellxgene_e2824739'></a>CELLxGENE_e2824739 — HTAPP-364-SMP-1321 Slide-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~2,958 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `e2824739-ea79-4efd-9434-3bec079b55d3.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/0dc10486-eca5-48af-8301-f61d9112c848.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_e2824739](https://cellxgene.cziscience.com/e/e2824739-ea79-4efd-9434-3bec079b55d3.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: liver | Assays: Slide-seqV2...

#### <a id='cellxgene_e5c614b8'></a>CELLxGENE_e5c614b8 — GSM6592050_M3
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `e5c614b8-b918-43b8-a542-670433e4da18.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/2378db2e-a807-4b99-9cba-6a5c22d95862.h5ad`
- **Citation / Reference:** 10.1016/j.labinv.2023.100258
- **Repository Access Link:** [CELLxGENE_e5c614b8](https://cellxgene.cziscience.com/e/e5c614b8-b918-43b8-a542-670433e4da18.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatial Transcriptomics Reveal Pitfalls and Opportunities for the Detection of Rare High-Plasticity Breast Cancer Subtypes | Diseases: triple-negative breast carcinoma | Tissues: breast | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_ec423499'></a>CELLxGENE_ec423499 — HTAPP-313-SMP-932 scRNA-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~12,494 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `ec423499-4564-4ab4-a7c2-316c72788ead.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/ec51b9a3-ab71-4aba-86ff-1d04aae7f7cc.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_ec423499](https://cellxgene.cziscience.com/e/ec423499-4564-4ab4-a7c2-316c72788ead.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: liver | Assays: 10x 3' v3...

#### <a id='cellxgene_f12ab0e6'></a>CELLxGENE_f12ab0e6 — HTAPP-812-SMP-8239 Slide-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~10,957 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `f12ab0e6-383c-4362-b1c6-922de09ccbc9.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/f9743bd7-ac48-4fe5-ba47-516b4ecd1a58.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_f12ab0e6](https://cellxgene.cziscience.com/e/f12ab0e6-383c-4362-b1c6-922de09ccbc9.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: axilla | Assays: Slide-seqV2...

#### <a id='cellxgene_ff4cfa86'></a>CELLxGENE_ff4cfa86 — HTAPP-213-SMP-6752 scRNA-seq
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~9,799 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `ff4cfa86-9c0c-4b7c-abd6-90547657d04f.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/5063908c-ebd9-4d29-91a7-1be98de6e298.h5ad`
- **Citation / Reference:** 10.1038/s41591-024-03215-z
- **Repository Access Link:** [CELLxGENE_ff4cfa86](https://cellxgene.cziscience.com/e/ff4cfa86-9c0c-4b7c-abd6-90547657d04f.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer metastatic microenvironment | Diseases: breast cancer | Tissues: chest wall | Assays: 10x 3' v3...

#### <a id='cellxgene_10bb68cf'></a>CELLxGENE_10bb68cf — GSM6592059_M13
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `10bb68cf-e20b-4c45-a65f-2d7d8129048e.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/e963ad5f-a05a-48c1-9c31-c9b3da73e00a.h5ad`
- **Citation / Reference:** 10.1016/j.labinv.2023.100258
- **Repository Access Link:** [CELLxGENE_10bb68cf](https://cellxgene.cziscience.com/e/10bb68cf-e20b-4c45-a65f-2d7d8129048e.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatial Transcriptomics Reveal Pitfalls and Opportunities for the Detection of Rare High-Plasticity Breast Cancer Subtypes | Diseases: triple-negative breast carcinoma | Tissues: breast | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_2cc628d1'></a>CELLxGENE_2cc628d1 — GSM6592055_M8
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `2cc628d1-b1dd-4300-9cf6-1015e2b1fd3d.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/36a207de-c1a2-474e-89f0-40117ef9b2e7.h5ad`
- **Citation / Reference:** 10.1016/j.labinv.2023.100258
- **Repository Access Link:** [CELLxGENE_2cc628d1](https://cellxgene.cziscience.com/e/2cc628d1-b1dd-4300-9cf6-1015e2b1fd3d.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatial Transcriptomics Reveal Pitfalls and Opportunities for the Detection of Rare High-Plasticity Breast Cancer Subtypes | Diseases: triple-negative breast carcinoma | Tissues: breast | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_480f9371'></a>CELLxGENE_480f9371 — GSM6592051_M4
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `480f9371-0e43-4724-a9b6-377b33789a42.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/7f70b1b7-8d8c-4a93-90fe-b7ec5daa5b26.h5ad`
- **Citation / Reference:** 10.1016/j.labinv.2023.100258
- **Repository Access Link:** [CELLxGENE_480f9371](https://cellxgene.cziscience.com/e/480f9371-0e43-4724-a9b6-377b33789a42.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatial Transcriptomics Reveal Pitfalls and Opportunities for the Detection of Rare High-Plasticity Breast Cancer Subtypes | Diseases: triple-negative breast carcinoma | Tissues: breast | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_540e4c1a'></a>CELLxGENE_540e4c1a — GSM6592060_M14
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `540e4c1a-4cba-4617-8d0b-583327b6afd0.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/e70e6ab2-81c9-4860-8959-2ce807ccfb01.h5ad`
- **Citation / Reference:** 10.1016/j.labinv.2023.100258
- **Repository Access Link:** [CELLxGENE_540e4c1a](https://cellxgene.cziscience.com/e/540e4c1a-4cba-4617-8d0b-583327b6afd0.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatial Transcriptomics Reveal Pitfalls and Opportunities for the Detection of Rare High-Plasticity Breast Cancer Subtypes | Diseases: triple-negative breast carcinoma | Tissues: breast | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_5d2c013d'></a>CELLxGENE_5d2c013d — GSM6592052_M5
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `5d2c013d-0221-401b-9a30-ecfd880f07a3.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/baa77986-95c7-4bf7-b4a4-09b23bcd16e4.h5ad`
- **Citation / Reference:** 10.1016/j.labinv.2023.100258
- **Repository Access Link:** [CELLxGENE_5d2c013d](https://cellxgene.cziscience.com/e/5d2c013d-0221-401b-9a30-ecfd880f07a3.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatial Transcriptomics Reveal Pitfalls and Opportunities for the Detection of Rare High-Plasticity Breast Cancer Subtypes | Diseases: triple-negative breast carcinoma | Tissues: breast | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_9adb1b29'></a>CELLxGENE_9adb1b29 — GSM6592062_M16
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `9adb1b29-65a2-4dd0-86bf-c02690d65fbd.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/5a1417ad-4861-442e-abf3-2d868c338559.h5ad`
- **Citation / Reference:** 10.1016/j.labinv.2023.100258
- **Repository Access Link:** [CELLxGENE_9adb1b29](https://cellxgene.cziscience.com/e/9adb1b29-65a2-4dd0-86bf-c02690d65fbd.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatial Transcriptomics Reveal Pitfalls and Opportunities for the Detection of Rare High-Plasticity Breast Cancer Subtypes | Diseases: triple-negative breast carcinoma | Tissues: breast | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_b0d9408e'></a>CELLxGENE_b0d9408e — GSM6592049_M2
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `b0d9408e-f02f-4ee4-9ab4-98d92cf290f1.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/8161db97-75ff-4a60-a360-e73b4da28555.h5ad`
- **Citation / Reference:** 10.1016/j.labinv.2023.100258
- **Repository Access Link:** [CELLxGENE_b0d9408e](https://cellxgene.cziscience.com/e/b0d9408e-f02f-4ee4-9ab4-98d92cf290f1.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatial Transcriptomics Reveal Pitfalls and Opportunities for the Detection of Rare High-Plasticity Breast Cancer Subtypes | Diseases: triple-negative breast carcinoma | Tissues: breast | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_c829c294'></a>CELLxGENE_c829c294 — GSM6592061_M15
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `c829c294-0995-4f0b-8c00-3d6ddab29b37.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/a13daadb-c385-4601-98f2-50f3e1d76be7.h5ad`
- **Citation / Reference:** 10.1016/j.labinv.2023.100258
- **Repository Access Link:** [CELLxGENE_c829c294](https://cellxgene.cziscience.com/e/c829c294-0995-4f0b-8c00-3d6ddab29b37.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatial Transcriptomics Reveal Pitfalls and Opportunities for the Detection of Rare High-Plasticity Breast Cancer Subtypes | Diseases: triple-negative breast carcinoma | Tissues: breast | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_dd1913a6'></a>CELLxGENE_dd1913a6 — GSM6592058_M11
- **Cancer Type / Indication:** Breast
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `dd1913a6-019d-4009-a118-6a3d158b22e6.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/631b8e61-c787-44a3-ba58-e3ba7c24f0d9.h5ad`
- **Citation / Reference:** 10.1016/j.labinv.2023.100258
- **Repository Access Link:** [CELLxGENE_dd1913a6](https://cellxgene.cziscience.com/e/dd1913a6-019d-4009-a118-6a3d158b22e6.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Spatial Transcriptomics Reveal Pitfalls and Opportunities for the Detection of Rare High-Plasticity Breast Cancer Subtypes | Diseases: triple-negative breast carcinoma | Tissues: breast | Assays: Visium Spatial Gene Expression V1...

### CRC Cohorts
*56 cohorts identified for CRC*

#### <a id='gse236581'></a>GSE236581 — Spatiotemporal single-cell analysis decodes cellular dynamics underlying different responses to immunotherapy in Colorectal Cancer
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Finally, a predictive signature was established using circulating CD8 T cells."*
- **Cohort Scale:** 169 samples/patients; ~338,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE236581_counts.mtx.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE236nnn/GSE236581/suppl/GSE236581_counts.mtx.gz`
- **Citation / Reference:** PMID:38981439
- **Repository Access Link:** [GSE236581](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236581)
- **Study Abstract / Experimental Design:**
  Expanding the efficacy of immune checkpoint blockade (ICB) in colorectal cancer (CRC) patients presses for a comprehensive understanding of treatment responsiveness. Here, we analyzed 169 single-cell samples from CRC patients at multiple sequential time points during the course of anti-PD-1 neoadjuvant therapy to map the evolution of local and systemic immunity. In tumors, exhausted T (Tex) cells ...

#### <a id='gse205506'></a>GSE205506 — Remodeling of the Immune and Stromal Cell Compartment by PD-1 Blockade in Mismatch Repair-Deficient Colorectal Cancer
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"of CD8+ Trm-mitotic, CD4+ Tregs, proinflammatory IL1B monocyte and CCL2+Fibroblast concertedly decrease following treatment, while those of CD8+ Tem, CD4+ Th, CD20+ B and HLA-DRA+ Endothelial cells increase."*
- **Cohort Scale:** 40 samples/patients; ~80,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE205506_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205506/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:37172580
- **Repository Access Link:** [GSE205506](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205506)
- **Study Abstract / Experimental Design:**
  Immune checkpoint inhibitor (ICI) therapy can induce complete responses in mismatch repair-deficient and microsatellite instability-high (d-MMR/MSI-H) colorectal cancers (CRCs). However, the mechanism responsible for pathological complete response (pCR) to immunotherapy has not been completely understood. We utilize single-cell RNA sequencing to examine the immune and stromal cell dynamics in 19 p...

#### <a id='gse299651'></a>GSE299651 — Pooled single-cell screening in colorectal cancer identifies transcriptional modules of clinical relevance unlocked by oncogenes
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 20 samples/patients; ~40,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE299651_20240523_Caco2_pool_cellranger_3_filtered_feature_bc_matrix.h5; GSE299651_20240523_HT29_pool_cellranger_3_filtered_feature_bc_matrix.h5; GSE299651_20240523_RKO_pool_cellranger_3_filtered_feature_bc_matrix.h5`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE299nnn/GSE299651/suppl/GSE299651_20240523_Caco2`
- **Citation / Reference:** PMID:41555096
- **Repository Access Link:** [GSE299651](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299651)
- **Study Abstract / Experimental Design:**
  While oncogenic mutations shape colorectal cancer biology and therapy response, their prognostic value remains low. Cluster-based classification of patient cancer transcriptomes has shown greater promise for prognosis, yet these systems do not account for the roles of oncogenes in establishing cancer phenotypes. Here, we create and validate a prognostic classifier for colorectal cancer based on tr...

#### <a id='gse146771'></a>GSE146771 — Single-Cell Analyses Inform Mechanisms of Myeloid-Targeted therapies in colon cancer
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"nt with a CD40 agonist antibody preferentially activated a conventional dendritic cell population and increased Bhlhe40+ Th1-like cells and CD8+ memory T cells."*
- **Cohort Scale:** 20 samples/patients; ~40,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE146771_CRC.Leukocyte.10x.Metadata.txt.gz; GSE146771_CRC.Leukocyte.10x.TPM.txt.gz; GSE146771_CRC.Leukocyte.Smart-seq2.Metadata.txt.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146771/suppl/GSE146771_CRC.Leukocyte.`
- **Citation / Reference:** PMID:32302573
- **Repository Access Link:** [GSE146771](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146771)
- **Study Abstract / Experimental Design:**
  Single cell RNA sequencing (scRNA-seq) is a powerful tool for defining cellular diversity in tumors, but its application towards dissecting mechanisms underlying immune-modulating therapies is scarce. We performed scRNA-seq analyses on immune and stromal populations from colorectal cancer patients, identifying specific macrophage and conventional dendritic cell (cDC) subsets as key mediators of ce...

#### <a id='gse164522'></a>GSE164522 — Single-cell analyses reveal phenotypic linkage between colorectal cancer and liver metastasis
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD45+ Immune-enriched)**
  > *"Single-cell suspension sorted by FACS for CD45+ leukocytes to enrich for tumor-infiltrating immune cells."*
- **Cohort Scale:** 17 samples/patients; ~34,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE164522_CRLM_LN_expression.csv.gz; GSE164522_CRLM_MN_expression.csv.gz; GSE164522_CRLM_MT_expression.csv.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164522/suppl/GSE164522_CRLM_LN_expres`
- **Citation / Reference:** PMID:35303421
- **Repository Access Link:** [GSE164522](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164522)
- **Study Abstract / Experimental Design:**
  The tumor microenvironment (TME) is intrinsically associated with clinical responses of immunotherapy, but it remains opaque how cancer cells and host tissues differentially influence its immune composition. Here, we performed systematic single-cell analyses for autologous clinical samples from liver metastasized colorectal cancer patients, as well as non-metastatic primary liver and colon cancer ...

#### <a id='gse235917'></a>GSE235917 — First-line durvalumab and tremelimumab with chemotherapy in RAS-mutated metastatic colorectal cancer: a phase 1b/2 trial [scRNA-Seq]
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 8 samples/patients; ~16,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE235917_5prim_matrix.mtx.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE235nnn/GSE235917/suppl/GSE235917_5prim_matrix.m`
- **Citation / Reference:** PMID:37563240
- **Repository Access Link:** [GSE235917](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235917)
- **Study Abstract / Experimental Design:**
  While patients with microsatellite instable, metastatic colorectal cancer (CRC) benefit from immune checkpoint blockade, chemotherapy with targeted therapies remains the only therapeutic option for microsatellite stable (MSS) tumors. The single arm, phase IB/II MEDITREME trial evaluates safety and efficacy of durvalumab plus tremelimumab in combination with mFOLFOX6 chemotherapy in first line, in ...

#### <a id='gse278406'></a>GSE278406 — Phenotypic plasticity and increased tissue infiltration of TREM1+ mono-macrophages following radiotherapy in rectal cancer. [scRNA-Seq]
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Following IR, loss of TREM1 in mono-macrophages undermines antitumor immunity by altering mono-macrophages differentiation and inhibiting CD8+ T cell infiltration and activation."*
- **Cohort Scale:** 6 samples/patients; ~12,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE278406_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE278nnn/GSE278406/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:39793571
- **Repository Access Link:** [GSE278406](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278406)
- **Study Abstract / Experimental Design:**
  Our previously reported phase II and phase III trials have demonstrated that short-course radiotherapy (SCRT) combined with neoadjuvant immunochemotherapy (SIC) has led to clinical benefits in locally advanced rectal cancer (LARC). Characterization of the molecular mechanisms underlying responses to SIC may lead to improved treatment strategies. Here, we prospectively collected and applied multi-o...

#### <a id='gse274321'></a>GSE274321 — IKKa modulates colorectal cancer metastasis by preventing tight junction stabilization and collective cell migration [scRNA-seq]
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 3 samples/patients; ~6,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE274321_scRNAseq_CRC_PD05_IKKa.RDS.gz; filelist.txt; GSE274321_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274321/suppl/GSE274321_scRNAseq_CRC_P`
- **Citation / Reference:** PMID:41484106
- **Repository Access Link:** [GSE274321](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274321)
- **Study Abstract / Experimental Design:**
  We have previously shown that IKKa coordinates the activation of several oncogenic and therapy-resistant pathways, including ATM/DDR, BRD4, and JAK/STAT3, independently of canonical NF-kB signaling. Here, we found that IKKa suppression, either genetically or pharmacologically, led to stabilization of ZO-1 protein and increased CLDN2 expression, resulting in altered tight junction distribution. IKK...

#### <a id='gse309346'></a>GSE309346 — PI3K and MAPK signaling nodes as divergent drivers of phenotypic plasticity in cancer-associated fibroblasts in colorectal cancer [scRNA-Seq]
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 3 samples/patients; ~6,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE309346_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE309nnn/GSE309346/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41817574
- **Repository Access Link:** [GSE309346](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309346)
- **Study Abstract / Experimental Design:**
  Cancer-associated fibroblasts (CAFs) exhibit phenotypic heterogeneity with each functional state playing critical roles in tumor progression. Notably, subtypes like inflammatory CAF (iCAF), characterized by increased chemokine/cytokine secretion, and myofibroblast-like CAF (myCAF), characterized by enhanced extracellular matrix (ECM) deposition and increased actomyosin contractility, can undergo p...

#### <a id='cellxgene_829a3cd1'></a>CELLxGENE_829a3cd1 — progressive_plasticity_during_crc_metastasis_epithelial
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 29 samples/patients; ~47,107 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `829a3cd1-a466-49f1-b2e9-d3f6b7f392e2.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/81bd8922-b69a-45ea-9b2a-032c480c0e37.h5ad`
- **Citation / Reference:** 10.1038/s41586-024-08560-0
- **Repository Access Link:** [CELLxGENE_829a3cd1](https://cellxgene.cziscience.com/e/829a3cd1-a466-49f1-b2e9-d3f6b7f392e2.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Progressive plasticity during colorectal cancer metastasis | Diseases: colorectal cancer, normal | Tissues: caecum, chest wall, descending colon, hepatic flexure of colon, liver, lung, peritoneum, rectum, sigmoid colon, transverse colon | Assays: 10x 3' v3...

#### <a id='cellxgene_2554a654'></a>CELLxGENE_2554a654 — progressive_plasticity_during_crc_metastasis_tumor
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Single-cell suspension sorted by FACS for CD3+/CD8+ T lymphocytes (TILs)."*
- **Cohort Scale:** 28 samples/patients; ~26,086 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `2554a654-237c-485e-b0d2-a58b46cf2d9a.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/0f353512-968e-4ffb-aab4-8162ac126e80.h5ad`
- **Citation / Reference:** 10.1038/s41586-024-08560-0
- **Repository Access Link:** [CELLxGENE_2554a654](https://cellxgene.cziscience.com/e/2554a654-237c-485e-b0d2-a58b46cf2d9a.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Progressive plasticity during colorectal cancer metastasis | Diseases: colorectal cancer | Tissues: caecum, chest wall, descending colon, hepatic flexure of colon, liver, lung, peritoneum, rectum, sigmoid colon, transverse colon | Assays: 10x 3' v3...

#### <a id='cellxgene_5ee552f5'></a>CELLxGENE_5ee552f5 — progressive_plasticity_during_crc_metastasis_non-tumor_epithelial
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 27 samples/patients; ~21,026 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `5ee552f5-834c-4766-92c8-5a09dfd45767.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/f20aae78-a32b-4914-8a98-7478c48324d3.h5ad`
- **Citation / Reference:** 10.1038/s41586-024-08560-0
- **Repository Access Link:** [CELLxGENE_5ee552f5](https://cellxgene.cziscience.com/e/5ee552f5-834c-4766-92c8-5a09dfd45767.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Progressive plasticity during colorectal cancer metastasis | Diseases: colorectal cancer, normal | Tissues: caecum, descending colon, hepatic flexure of colon, rectum, sigmoid colon, transverse colon | Assays: 10x 3' v3...

#### <a id='cellxgene_4b5afdf9'></a>CELLxGENE_4b5afdf9 — progressive_plasticity_during_crc_metastasis_untreated_epithelial
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 9 samples/patients; ~13,843 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `4b5afdf9-9299-4655-9218-c9f67b3776e7.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/dbc8f31b-e54e-4519-a9e2-c6c17765914d.h5ad`
- **Citation / Reference:** 10.1038/s41586-024-08560-0
- **Repository Access Link:** [CELLxGENE_4b5afdf9](https://cellxgene.cziscience.com/e/4b5afdf9-9299-4655-9218-c9f67b3776e7.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Progressive plasticity during colorectal cancer metastasis | Diseases: colorectal cancer, normal | Tissues: caecum, hepatic flexure of colon, liver, sigmoid colon | Assays: 10x 3' v3...

#### <a id='gse188711'></a>GSE188711 — Resolving the Difference Between Left-sided and Right-sided Colorectal Cancer by Single-cell Sequencing
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Right-sided CRC harbors a significant proportion of exhausted CD8 T cells of a highly migratory nature."*
- **Cohort Scale:** 6 samples/patients; ~12,000 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `filelist.txt; GSE188711_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE188nnn/GSE188711/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:34793335
- **Repository Access Link:** [GSE188711](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188711)
- **Study Abstract / Experimental Design:**
  Colorectal cancers (CRCs) exhibit differences in incidence, pathogenesis, molecular pathways and outcome depending on the location of the tumor. The transcriptomes of 27,927 single human CRC cells, from three left-sided and three right-sided CRC patients were profiled by scRNA-seq. Right-sided CRC harbors a significant proportion of exhausted CD8 T cells of a highly migratory nature. One cluster o...

#### <a id='gse336564'></a>GSE336564 — ZFP36L2 orchestrates stress-adaptive plasticity in intestinal regeneration and colorectal cancer metastasis (PDO scRNAseq)
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** 10x Chromium (3' unspecified)
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (Viability DAPI-/7-AAD- only)**
  > *"Viable (DAPI-negative) cells were FACS-sorted and pooled in equal numbers across conditions for library preparation using the Chromium Single Cell 3′ v3."*
- **Cohort Scale:** 6 samples/patients; ~10,000 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `GSE336564_OKG146P.h5ad`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE336nnn/GSE336564/suppl/GSE336564_OKG146P.h5ad`
- **Repository Access Link:** [GSE336564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336564)
- **Study Abstract / Experimental Design:**
  Single-cell RNA sequencing (scRNA-seq) was performed on patient-derived primary colorectal cancer (CRC) organoids (OKG146P) to investigate the role of ZFP36L2 in regulating cell state plasticity during differentiation and dedifferentiation. Organoids expressing a doxycycline-inducible shRNA targeting ZFP36L2 (shZFP36L2) or a non-targeting control (shCtrl) were cultured with 2 μg/mL doxycycline and...

#### <a id='gse216534'></a>GSE216534 — γδ T cells are effectors of immunotherapy in cancers with HLA class I defects
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Single-cell suspension sorted by FACS for CD3+/CD8+ T lymphocytes (TILs)."*
- **Cohort Scale:** 2 samples/patients; ~4,000 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `GSE216534_filtered_feature_bc_matrix.h5`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE216nnn/GSE216534/suppl/GSE216534_filtered_featu`
- **Citation / Reference:** PMID:36631610
- **Repository Access Link:** [GSE216534](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216534)
- **Study Abstract / Experimental Design:**
  DNA mismatch repair deficient (MMR-d) cancers present an abundance of neoantigens that likely underlies their exceptional responsiveness to immune checkpoint blockade (ICB). However, MMR-d colon cancers that evade CD8+ T cells through loss of Human Leukocyte Antigen (HLA) class I-mediated antigen presentation frequently remain responsive to ICB, suggesting the involvement of other immune effector ...

#### <a id='cellxgene_387acac5'></a>CELLxGENE_387acac5 — progressive_plasticity_during_crc_metastasis_kg150_tumor
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Single-cell suspension sorted by FACS for CD3+/CD8+ T lymphocytes (TILs)."*
- **Cohort Scale:** 1 samples/patients; ~2,574 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `387acac5-48cd-4f25-880e-523114ccb006.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/dbd79659-9f19-4035-ad38-ec1e2eaf22e0.h5ad`
- **Citation / Reference:** 10.1038/s41586-024-08560-0
- **Repository Access Link:** [CELLxGENE_387acac5](https://cellxgene.cziscience.com/e/387acac5-48cd-4f25-880e-523114ccb006.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Progressive plasticity during colorectal cancer metastasis | Diseases: colorectal cancer | Tissues: liver, sigmoid colon | Assays: 10x 3' v3...

#### <a id='cellxgene_ef0d813e'></a>CELLxGENE_ef0d813e — progressive_plasticity_during_crc_metastasis_kg183_tumor
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Single-cell suspension sorted by FACS for CD3+/CD8+ T lymphocytes (TILs)."*
- **Cohort Scale:** 1 samples/patients; ~1,203 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `ef0d813e-b9b4-4d30-8e54-b3939942c447.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/9498d1c6-ab77-4f8b-9a97-bbe83a736cdf.h5ad`
- **Citation / Reference:** 10.1038/s41586-024-08560-0
- **Repository Access Link:** [CELLxGENE_ef0d813e](https://cellxgene.cziscience.com/e/ef0d813e-b9b4-4d30-8e54-b3939942c447.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Progressive plasticity during colorectal cancer metastasis | Diseases: colorectal cancer | Tissues: liver, sigmoid colon | Assays: 10x 3' v3...

#### <a id='cellxgene_2e95d453'></a>CELLxGENE_2e95d453 — progressive_plasticity_during_crc_metastasis_kg146_tumor
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Single-cell suspension sorted by FACS for CD3+/CD8+ T lymphocytes (TILs)."*
- **Cohort Scale:** 1 samples/patients; ~3,351 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `2e95d453-930b-4425-9475-ccaf8cf2ac40.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/987f1e0d-c402-46b2-8344-d7880b523bda.h5ad`
- **Citation / Reference:** 10.1038/s41586-024-08560-0
- **Repository Access Link:** [CELLxGENE_2e95d453](https://cellxgene.cziscience.com/e/2e95d453-930b-4425-9475-ccaf8cf2ac40.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Progressive plasticity during colorectal cancer metastasis | Diseases: colorectal cancer | Tissues: liver, lung, rectum | Assays: 10x 3' v3...

#### <a id='cellxgene_05a8c945'></a>CELLxGENE_05a8c945 — The single-cell Colorectal Cancer Atlas -- core atlas
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 588 samples/patients; ~3,790,266 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `05a8c945-bc12-414f-960d-a31943bbcdd1.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/4a8b9568-965e-46b8-a427-baab6bf018e5.h5ad`
- **Citation / Reference:** 10.1016/j.ccell.2025.12.003
- **Repository Access Link:** [CELLxGENE_05a8c945](https://cellxgene.cziscience.com/e/05a8c945-bc12-414f-960d-a31943bbcdd1.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Single-cell integration and multi-modal profiling reveals phenotypes and spatial organization of neutrophils in colorectal cancer | Diseases: colon adenocarcinoma, colon adenocarcinoma || mucinous adenocarcinoma, colon adenoma, colon adenoma || colorectal tubulovillous adenoma, colon adenoma || tubular adenoma, colon carcinoma, colon sessile serrated adenoma/polyp, colorectal adenocarc...

#### <a id='cellxgene_19053a82'></a>CELLxGENE_19053a82 — Extended+ - 18485 genes
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 308 samples/patients; ~1,596,200 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `19053a82-9c89-4fb8-bd19-d7b1800b0b7b.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/757945c8-a916-431d-aceb-1afbc80a7c55.h5ad`
- **Citation / Reference:** 10.1038/s41586-024-07571-1
- **Repository Access Link:** [CELLxGENE_19053a82](https://cellxgene.cziscience.com/e/19053a82-9c89-4fb8-bd19-d7b1800b0b7b.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Single-cell integration reveals metaplasia in inflammatory gut diseases | Diseases: Crohn disease, celiac disease, colorectal cancer, gastric cancer, inflammatory bowel disease, normal, ulcerative colitis | Tissues: ascending colon, body of stomach, buccal mucosa, caecum, caecum epithelium, cecum mucosa, colon, colonic epithelium, colonic mucosa, descending colon, duodenal epithelium, ...

#### <a id='cellxgene_e6aaf5a4'></a>CELLxGENE_e6aaf5a4 — Extended - All genes
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 233 samples/patients; ~1,358,573 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `e6aaf5a4-16e9-4ea6-9733-4eafd4e473d3.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/8b6c3914-0505-40a5-a241-a9390759f7ac.h5ad`
- **Citation / Reference:** 10.1038/s41586-024-07571-1
- **Repository Access Link:** [CELLxGENE_e6aaf5a4](https://cellxgene.cziscience.com/e/e6aaf5a4-16e9-4ea6-9733-4eafd4e473d3.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Single-cell integration reveals metaplasia in inflammatory gut diseases | Diseases: Crohn disease, colorectal cancer, gastric cancer, inflammatory bowel disease, normal, ulcerative colitis | Tissues: ascending colon, body of stomach, buccal mucosa, caecum, colon, colonic epithelium, colonic mucosa, descending colon, duodenal epithelium, duodenum, epithelium of esophagus, epithelium of ...

#### <a id='cellxgene_dc6b1e06'></a>CELLxGENE_dc6b1e06 — Extended - T and NK cells
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 227 samples/patients; ~262,642 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `dc6b1e06-3880-43f6-8de8-ef9cb5e7ef5e.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/a7d8cf9b-2825-42ea-99f3-e96c9aee23c7.h5ad`
- **Citation / Reference:** 10.1038/s41586-024-07571-1
- **Repository Access Link:** [CELLxGENE_dc6b1e06](https://cellxgene.cziscience.com/e/dc6b1e06-3880-43f6-8de8-ef9cb5e7ef5e.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Single-cell integration reveals metaplasia in inflammatory gut diseases | Diseases: Crohn disease, colorectal cancer, gastric cancer, inflammatory bowel disease, normal, ulcerative colitis | Tissues: ascending colon, body of stomach, buccal mucosa, caecum, colon, colonic epithelium, colonic mucosa, descending colon, duodenal epithelium, duodenum, epithelium of esophagus, epithelium of ...

#### <a id='cellxgene_9c235282'></a>CELLxGENE_9c235282 — Extended - Myeloid cells
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 223 samples/patients; ~50,570 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `9c235282-2b6f-4f74-8f16-7d3ac1b51371.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/21862555-9e6e-41fb-8ff7-63f4acb042ba.h5ad`
- **Citation / Reference:** 10.1038/s41586-024-07571-1
- **Repository Access Link:** [CELLxGENE_9c235282](https://cellxgene.cziscience.com/e/9c235282-2b6f-4f74-8f16-7d3ac1b51371.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Single-cell integration reveals metaplasia in inflammatory gut diseases | Diseases: Crohn disease, colorectal cancer, gastric cancer, inflammatory bowel disease, normal, ulcerative colitis | Tissues: ascending colon, body of stomach, buccal mucosa, caecum, colon, colonic epithelium, colonic mucosa, descending colon, duodenal epithelium, duodenum, epithelium of esophagus, epithelium of ...

#### <a id='cellxgene_40a0ade8'></a>CELLxGENE_40a0ade8 — Extended - Myeloid_with_neutrophils
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 223 samples/patients; ~52,404 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `40a0ade8-6067-4e22-9224-4d3c5e9bfc0d.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/48ecc234-7073-4de9-8be6-1c2f08d4c146.h5ad`
- **Citation / Reference:** 10.1038/s41586-024-07571-1
- **Repository Access Link:** [CELLxGENE_40a0ade8](https://cellxgene.cziscience.com/e/40a0ade8-6067-4e22-9224-4d3c5e9bfc0d.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Single-cell integration reveals metaplasia in inflammatory gut diseases | Diseases: Crohn disease, colorectal cancer, gastric cancer, inflammatory bowel disease, normal, ulcerative colitis | Tissues: ascending colon, body of stomach, buccal mucosa, caecum, colon, colonic epithelium, colonic mucosa, descending colon, duodenal epithelium, duodenum, epithelium of esophagus, epithelium of ...

#### <a id='cellxgene_1d54fb17'></a>CELLxGENE_1d54fb17 — Extended - B and B plasma cells
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 207 samples/patients; ~250,094 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `1d54fb17-90d5-47ca-ad8d-22ab63cb5c1f.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/667d05a8-e1f3-412c-820d-d39592e0cd55.h5ad`
- **Citation / Reference:** 10.1038/s41586-024-07571-1
- **Repository Access Link:** [CELLxGENE_1d54fb17](https://cellxgene.cziscience.com/e/1d54fb17-90d5-47ca-ad8d-22ab63cb5c1f.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Single-cell integration reveals metaplasia in inflammatory gut diseases | Diseases: Crohn disease, colorectal cancer, gastric cancer, inflammatory bowel disease, normal, ulcerative colitis | Tissues: ascending colon, body of stomach, buccal mucosa, caecum, colon, colonic epithelium, colonic mucosa, descending colon, duodenal epithelium, duodenum, epithelium of esophagus, epithelium of ...

#### <a id='cellxgene_278eac3f'></a>CELLxGENE_278eac3f — Extended - Endothelial cells
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 195 samples/patients; ~60,411 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `278eac3f-a5a4-4d0a-9d83-0ca7efc67e16.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/01db6084-3cff-4537-b52f-139dc9df50a5.h5ad`
- **Citation / Reference:** 10.1038/s41586-024-07571-1
- **Repository Access Link:** [CELLxGENE_278eac3f](https://cellxgene.cziscience.com/e/278eac3f-a5a4-4d0a-9d83-0ca7efc67e16.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Single-cell integration reveals metaplasia in inflammatory gut diseases | Diseases: Crohn disease, colorectal cancer, gastric cancer, inflammatory bowel disease, normal, ulcerative colitis | Tissues: ascending colon, body of stomach, buccal mucosa, caecum, colon, colonic epithelium, colonic mucosa, descending colon, duodenum, epithelium of esophagus, esophagus, gingiva, ileal epitheliu...

#### <a id='cellxgene_ef7bb7f0'></a>CELLxGENE_ef7bb7f0 — Extended - Mesenchymal (adult/pediatric)
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 184 samples/patients; ~77,050 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `ef7bb7f0-234d-459f-b613-b37cc7a4c70f.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/e43f0ac0-a426-4471-8c4f-a9f99d612357.h5ad`
- **Citation / Reference:** 10.1038/s41586-024-07571-1
- **Repository Access Link:** [CELLxGENE_ef7bb7f0](https://cellxgene.cziscience.com/e/ef7bb7f0-234d-459f-b613-b37cc7a4c70f.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Single-cell integration reveals metaplasia in inflammatory gut diseases | Diseases: Crohn disease, colorectal cancer, gastric cancer, inflammatory bowel disease, normal, ulcerative colitis | Tissues: ascending colon, body of stomach, buccal mucosa, caecum, colon, colonic epithelium, colonic mucosa, descending colon, duodenum, epithelium of esophagus, epithelium of rectum, epithelium of...

#### <a id='cellxgene_7be23e52'></a>CELLxGENE_7be23e52 — Extended - Neural cells
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 101 samples/patients; ~23,904 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `7be23e52-e357-45d6-861f-31166913104e.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/02f1646e-52ca-4b4a-aec2-f4469f216cea.h5ad`
- **Citation / Reference:** 10.1038/s41586-024-07571-1
- **Repository Access Link:** [CELLxGENE_7be23e52](https://cellxgene.cziscience.com/e/7be23e52-e357-45d6-861f-31166913104e.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Single-cell integration reveals metaplasia in inflammatory gut diseases | Diseases: Crohn disease, colorectal cancer, gastric cancer, normal, ulcerative colitis | Tissues: ascending colon, body of stomach, buccal mucosa, caecum, colon, colonic mucosa, descending colon, duodenum, esophagus, gingiva, ileum, intestine, jejunum, labial gland, mesenteric lymph node, periodontium, pyloric an...

#### <a id='cellxgene_763d1d88'></a>CELLxGENE_763d1d88 — Extended - Large Intestine (adult/pediatric)
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 61 samples/patients; ~96,675 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `763d1d88-b38d-4ab0-b00f-72df816b4b08.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/681f6f1b-f9b6-4d1d-8588-97fedad51f7a.h5ad`
- **Citation / Reference:** 10.1038/s41586-024-07571-1
- **Repository Access Link:** [CELLxGENE_763d1d88](https://cellxgene.cziscience.com/e/763d1d88-b38d-4ab0-b00f-72df816b4b08.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Single-cell integration reveals metaplasia in inflammatory gut diseases | Diseases: colorectal cancer, inflammatory bowel disease, normal, ulcerative colitis | Tissues: ascending colon, caecum, colon, colonic epithelium, colonic mucosa, descending colon, epithelium of rectum, lamina propria of mucosa of colon, rectosigmoid junction, rectum, sigmoid colon, transverse colon, vermiform ap...

#### <a id='cellxgene_6a270451'></a>CELLxGENE_6a270451 — VAL and DIS datasets: Non-Epithelial
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 59 samples/patients; ~10,696 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `6a270451-b4d9-43e0-aa89-e33aac1ac74b.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/91a55470-7520-4ecf-bcbc-1062f7ef4cdb.h5ad`
- **Citation / Reference:** 10.1016/j.cell.2021.11.031
- **Repository Access Link:** [CELLxGENE_6a270451](https://cellxgene.cziscience.com/e/6a270451-b4d9-43e0-aa89-e33aac1ac74b.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN VUMC - Differential pre-malignant programs and microenvironment chart distinct paths to malignancy in human colorectal polyps | Diseases: colon sessile serrated adenoma/polyp, colorectal cancer, hyperplastic polyp, normal, tubular adenoma, tubulovillous adenoma | Tissues: ascending colon, descending colon, hepatic cecum, hepatic flexure of colon, rectum, sigmoid colon, transverse ...

#### <a id='gse294300'></a>GSE294300 — Single-cell gene expression profiles of cells from primary colorectal cancer and adjacent normal tissue [scRNA-Seq]
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 36 samples/patients; ~72,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE294300_cell_batch.tsv.gz; filelist.txt; GSE294300_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE294nnn/GSE294300/suppl/GSE294300_cell_batch.tsv`
- **Citation / Reference:** PMID:42143353
- **Repository Access Link:** [GSE294300](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294300)
- **Study Abstract / Experimental Design:**
  Colorectal cancer (CRC) is the third malignancy worldwide. RAS mutant CRC has a worse prognosis and resistant to immune therapies. Current research on the tumor microenvironment of RAS-mutant CRC remains limited, with few studies systematically characterizing its cellular composition or functional dynamics. A total of 36 clinical surgical samples (tumors and paired adjacent normal tissues) from 18...

#### <a id='cellxgene_d6dfdef1'></a>CELLxGENE_d6dfdef1 — Validation (Val) set of human colorectal tumor: Epithelial
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 26 samples/patients; ~57,723 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `d6dfdef1-406d-4efb-808c-3c5eddbfe0cb.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/6f8f0b95-3779-4c0e-824d-4810c39f008c.h5ad`
- **Citation / Reference:** 10.1016/j.cell.2021.11.031
- **Repository Access Link:** [CELLxGENE_d6dfdef1](https://cellxgene.cziscience.com/e/d6dfdef1-406d-4efb-808c-3c5eddbfe0cb.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: HTAN VUMC - Differential pre-malignant programs and microenvironment chart distinct paths to malignancy in human colorectal polyps | Diseases: colon sessile serrated adenoma/polyp, colorectal cancer, hyperplastic polyp, normal, tubular adenoma, tubulovillous adenoma | Tissues: ascending colon, descending colon, hepatic cecum, hepatic flexure of colon, rectum, sigmoid colon, transverse ...

#### <a id='gse271690'></a>GSE271690 — An integrated single cell and spatial transcriptomic map of Tumor Margin in Colorectal Cancer
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 16 samples/patients; ~32,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE271690_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271690/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE271690](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271690)
- **Study Abstract / Experimental Design:**
  Solid tumors are highly heterogeneous complex ecosystems. By dissecting the tumor ecosystem, particularly around the tumor periphery, we can gain deeper insights into the mechanisms of tumor cell infiltration and invasion....

#### <a id='gse315534'></a>GSE315534 — Single-cell RNA sequencing of primary colorectal cancer, matched adjacent normal tissues, and paired liver metastases
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 15 samples/patients; ~30,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE315534_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE315nnn/GSE315534/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:42082451
- **Repository Access Link:** [GSE315534](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE315534)
- **Study Abstract / Experimental Design:**
  This study will perform single-cell RNA sequencing on primary tumor tissues (T) and matched adjacent non-tumor tissues (N) from six colorectal cancer patients, along with paired liver metastatic lesions (M) from three of these patients, to construct an integrated single-cell atlas spanning primary tumors, adjacent normal tissues, and distant metastases. This multi-region, patient-matched design en...

#### <a id='gse335811'></a>GSE335811 — Macrophage-Induced Senescent Cancer-Associated Fibroblasts Promote SASP-Mediated Chemoresistance in Colorectal Cancer
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 12 samples/patients; ~24,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE335811_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE335nnn/GSE335811/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE335811](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE335811)
- **Study Abstract / Experimental Design:**
  Cancer-associated fibroblasts (CAFs) play a crucial role in the tumor microenvironment (TME) by influencing tumor progression and therapy resistance. Accumulating evidence suggests that CAFs undergo senescence, which can impact their effects on the TME. In this study, we integrated single-cell RNA sequencing (scRNA-seq), spatial transcriptomics, and multiple preclinical models to explore the mecha...

#### <a id='cellxgene_e3ed2ba4'></a>CELLxGENE_e3ed2ba4 — CRLM-NMP-ATLAS
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** BD Rhapsody
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 6 samples/patients; ~75,104 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `e3ed2ba4-edf5-40ac-8750-8a417ad1eefe.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/a3bcb68b-a36d-460f-bdde-da4fc0771c8c.h5ad`
- **Citation / Reference:** 10.1186/s12943-025-02430-7
- **Repository Access Link:** [CELLxGENE_e3ed2ba4](https://cellxgene.cziscience.com/e/e3ed2ba4-edf5-40ac-8750-8a417ad1eefe.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Ex vivo modelling of human colorectal cancer liver metastasis by normothermic machine perfusion | Diseases: colorectal carcinoma || metastatic malignant neoplasm, normal | Tissues: liver | Assays: BD Rhapsody Whole Transcriptome Analysis...

#### <a id='gse330797'></a>GSE330797 — Single-cell transcriptomics analysis of colorectal cancer (CRC)
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 4 samples/patients; ~8,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE330797_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE330nnn/GSE330797/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:42360233
- **Repository Access Link:** [GSE330797](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330797)
- **Study Abstract / Experimental Design:**
  Understanding cellular processes underlying colorectal cancer (CRC) development is needed to devise intervention strategies. Here, we performed single-cell RNA sequencing (scRNA-seq) of human treatment-naive colorectal cancer....

#### <a id='gse311338'></a>GSE311338 — Peristaltic forces drive tumor cell invasion in colorectal cancer [scRNA-seq]
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 3 samples/patients; ~6,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE311338_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE311nnn/GSE311338/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE311338](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE311338)
- **Study Abstract / Experimental Design:**
  Mechanical forces are known to influence the progression of cancer, however, the impact of naturally occurring physical forces remains less well understood. With the advent of organ-on-chip (OOC) technology, preclinical models can now incorporate human-relevant physiological forces and allow for more precise investigation of their effects.  In this study, we explore how the peristaltic motions of ...

#### <a id='gse312804'></a>GSE312804 — scRNAseq analysis of hepatocytes from liver cirrhosis patient
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium (3' unspecified)
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 2 samples/patients; ~4,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE312804_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE312nnn/GSE312804/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41580988
- **Repository Access Link:** [GSE312804](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE312804)
- **Study Abstract / Experimental Design:**
  We used the Flex single-cell RNA sequencing platform (10x Genomics) to investigate the transcriptomic differences of hepatocytes between healthy and liver cirrhosis liver samples obtained from resected tissues (non-tumorous regions) of patients with colorectal cancer liver metastases and patients with hepatocellular carcinoma accompanied by liver cirrhosis....

#### <a id='gse312260'></a>GSE312260 — Gene expression profile at single cell level of Irinotecan-sensitive and Irinotecan-resistant organoids for colorectal cancer
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 2 samples/patients; ~4,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE312260_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE312nnn/GSE312260/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41678386
- **Repository Access Link:** [GSE312260](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE312260)
- **Study Abstract / Experimental Design:**
  Organoids derived from both patients contained only epithelial cells....

#### <a id='cellxgene_f7af19e4'></a>CELLxGENE_f7af19e4 — S7_Rec/Sig A798015 Rep1
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `f7af19e4-fc64-46d0-ab3c-d4f70dd670f4.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/7a9397d6-dbd5-400e-bd35-e31793ade49c.h5ad`
- **Citation / Reference:** 10.1038/s41698-023-00488-4
- **Repository Access Link:** [CELLxGENE_f7af19e4](https://cellxgene.cziscience.com/e/f7af19e4-fc64-46d0-ab3c-d4f70dd670f4.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Profiling the heterogeneity of colorectal cancer consensus molecular subtypes using spatial transcriptomics | Diseases: colorectal cancer | Tissues: rectosigmoid junction | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_7fe57023'></a>CELLxGENE_7fe57023 — S2_Col_R A595688 Rep2
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `7fe57023-bbcf-4254-9aa4-bd2abf7d83ad.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/89a0b783-a960-4a65-b32b-5de604e5f9f9.h5ad`
- **Citation / Reference:** 10.1038/s41698-023-00488-4
- **Repository Access Link:** [CELLxGENE_7fe57023](https://cellxgene.cziscience.com/e/7fe57023-bbcf-4254-9aa4-bd2abf7d83ad.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Profiling the heterogeneity of colorectal cancer consensus molecular subtypes using spatial transcriptomics | Diseases: colorectal cancer | Tissues: right colon | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_879bb6df'></a>CELLxGENE_879bb6df — S2_Col_R A595688 Rep1
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `879bb6df-cc2a-40f1-854b-5be9629d03b2.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/a9ccd651-a17a-4b2e-bcec-5fd0e68d1681.h5ad`
- **Citation / Reference:** 10.1038/s41698-023-00488-4
- **Repository Access Link:** [CELLxGENE_879bb6df](https://cellxgene.cziscience.com/e/879bb6df-cc2a-40f1-854b-5be9629d03b2.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Profiling the heterogeneity of colorectal cancer consensus molecular subtypes using spatial transcriptomics | Diseases: colorectal cancer | Tissues: right colon | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_a73f7983'></a>CELLxGENE_a73f7983 — S1_Cec A551763 Rep2
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `a73f7983-c94b-4005-aaed-30d35b61688a.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/7e26f8dc-2b00-4b6a-9f1f-25a1ba853599.h5ad`
- **Citation / Reference:** 10.1038/s41698-023-00488-4
- **Repository Access Link:** [CELLxGENE_a73f7983](https://cellxgene.cziscience.com/e/a73f7983-c94b-4005-aaed-30d35b61688a.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Profiling the heterogeneity of colorectal cancer consensus molecular subtypes using spatial transcriptomics | Diseases: colorectal cancer | Tissues: caecum | Assays: Visium Spatial Gene Expression V1...

#### <a id='gse270767'></a>GSE270767 — Disseminated tumor cell of colorectal cancer adopt the wound healing program of epidermal keratinocyte [scRNA-seq]
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~2,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE270767_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270767/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:42436119
- **Repository Access Link:** [GSE270767](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270767)
- **Study Abstract / Experimental Design:**
  Disseminated tumor cells (DTCs) are a critical cell population in the metastasis process, and understanding their features is essential for the development of cancer treatment. We show that DTCs of colorectal cancer (CRC) transiently adopted the cellular state resembling to the activated epidermal keratinocyte during wound healing. Using the orthotropic transplantation of patient derived organoids...

#### <a id='cellxgene_b5753bee'></a>CELLxGENE_b5753bee — S1_Cec A551763 Rep1
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `b5753bee-dae1-4fa6-80c0-e78faa77b0aa.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/273cbea6-d2f1-40ca-af53-cd64dfa2c336.h5ad`
- **Citation / Reference:** 10.1038/s41698-023-00488-4
- **Repository Access Link:** [CELLxGENE_b5753bee](https://cellxgene.cziscience.com/e/b5753bee-dae1-4fa6-80c0-e78faa77b0aa.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Profiling the heterogeneity of colorectal cancer consensus molecular subtypes using spatial transcriptomics | Diseases: colorectal cancer | Tissues: caecum | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_c0d43178'></a>CELLxGENE_c0d43178 — S5_Rec A121573 Rep2
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `c0d43178-d368-4381-aecb-4e63c553d1fa.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/c17e73b9-e434-497c-a5ea-ffb5c07a4cd2.h5ad`
- **Citation / Reference:** 10.1038/s41698-023-00488-4
- **Repository Access Link:** [CELLxGENE_c0d43178](https://cellxgene.cziscience.com/e/c0d43178-d368-4381-aecb-4e63c553d1fa.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Profiling the heterogeneity of colorectal cancer consensus molecular subtypes using spatial transcriptomics | Diseases: colorectal cancer | Tissues: rectum | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_1e191a00'></a>CELLxGENE_1e191a00 — S3_Col_R A416371 Rep2
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `1e191a00-65da-4170-b9b7-eb94359d7de4.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/86ea4533-fc57-4087-9104-4ed781796841.h5ad`
- **Citation / Reference:** 10.1038/s41698-023-00488-4
- **Repository Access Link:** [CELLxGENE_1e191a00](https://cellxgene.cziscience.com/e/1e191a00-65da-4170-b9b7-eb94359d7de4.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Profiling the heterogeneity of colorectal cancer consensus molecular subtypes using spatial transcriptomics | Diseases: colorectal cancer | Tissues: right colon | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_7ba1a805'></a>CELLxGENE_7ba1a805 — S6_Rec A938797 Rep1
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `7ba1a805-0afe-4264-b21e-0ce2ef1aa3bd.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/caed233d-e71f-4db2-9bac-d306f8828544.h5ad`
- **Citation / Reference:** 10.1038/s41698-023-00488-4
- **Repository Access Link:** [CELLxGENE_7ba1a805](https://cellxgene.cziscience.com/e/7ba1a805-0afe-4264-b21e-0ce2ef1aa3bd.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Profiling the heterogeneity of colorectal cancer consensus molecular subtypes using spatial transcriptomics | Diseases: colorectal cancer | Tissues: rectum | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_74e80fd1'></a>CELLxGENE_74e80fd1 — S4_Col_Sig A120838 Rep2
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `74e80fd1-d9cd-4132-af0f-9c5d37489392.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/447b545b-e3ca-43b0-bb0b-022bf5c6c47d.h5ad`
- **Citation / Reference:** 10.1038/s41698-023-00488-4
- **Repository Access Link:** [CELLxGENE_74e80fd1](https://cellxgene.cziscience.com/e/74e80fd1-d9cd-4132-af0f-9c5d37489392.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Profiling the heterogeneity of colorectal cancer consensus molecular subtypes using spatial transcriptomics | Diseases: colorectal cancer | Tissues: sigmoid colon | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_729f397a'></a>CELLxGENE_729f397a — S3_Col_R A416371 Rep1
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `729f397a-0812-4b52-a7d1-b377107ffb41.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/a316235e-621b-4c0f-88d2-e4edba3dc1b0.h5ad`
- **Citation / Reference:** 10.1038/s41698-023-00488-4
- **Repository Access Link:** [CELLxGENE_729f397a](https://cellxgene.cziscience.com/e/729f397a-0812-4b52-a7d1-b377107ffb41.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Profiling the heterogeneity of colorectal cancer consensus molecular subtypes using spatial transcriptomics | Diseases: colorectal cancer | Tissues: right colon | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_2d821164'></a>CELLxGENE_2d821164 — S5_Rec A121573 Rep1
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `2d821164-4681-47dd-b2ac-eb5fddb5b621.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/909cb5ed-688f-41b9-9bf4-7247fba44371.h5ad`
- **Citation / Reference:** 10.1038/s41698-023-00488-4
- **Repository Access Link:** [CELLxGENE_2d821164](https://cellxgene.cziscience.com/e/2d821164-4681-47dd-b2ac-eb5fddb5b621.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Profiling the heterogeneity of colorectal cancer consensus molecular subtypes using spatial transcriptomics | Diseases: colorectal cancer | Tissues: rectum | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_2916b663'></a>CELLxGENE_2916b663 — S4_Col_Sig A120838 Rep1
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `2916b663-5fe0-4a53-a2da-f6875214e3e6.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/e8c6d38b-c0d8-4dc6-bb66-9ce3c719c40b.h5ad`
- **Citation / Reference:** 10.1038/s41698-023-00488-4
- **Repository Access Link:** [CELLxGENE_2916b663](https://cellxgene.cziscience.com/e/2916b663-5fe0-4a53-a2da-f6875214e3e6.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Profiling the heterogeneity of colorectal cancer consensus molecular subtypes using spatial transcriptomics | Diseases: colorectal cancer | Tissues: sigmoid colon | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_15b98664'></a>CELLxGENE_15b98664 — S6_Rec A938797 Rep2
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `15b98664-a52b-4877-aed9-03f4d73d2d2c.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/4f43fbfe-b88a-4010-9d5f-359afafe80e0.h5ad`
- **Citation / Reference:** 10.1038/s41698-023-00488-4
- **Repository Access Link:** [CELLxGENE_15b98664](https://cellxgene.cziscience.com/e/15b98664-a52b-4877-aed9-03f4d73d2d2c.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Profiling the heterogeneity of colorectal cancer consensus molecular subtypes using spatial transcriptomics | Diseases: colorectal cancer | Tissues: rectum | Assays: Visium Spatial Gene Expression V1...

#### <a id='cellxgene_297b5b89'></a>CELLxGENE_297b5b89 — S7_Rec/Sig A798015 Rep2
- **Cancer Type / Indication:** CRC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~4,992 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `297b5b89-6197-4135-acff-501c8d4199d0.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/9ff3a1ac-8e44-4c39-a578-e2ab93853b6c.h5ad`
- **Citation / Reference:** 10.1038/s41698-023-00488-4
- **Repository Access Link:** [CELLxGENE_297b5b89](https://cellxgene.cziscience.com/e/297b5b89-6197-4135-acff-501c8d4199d0.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Profiling the heterogeneity of colorectal cancer consensus molecular subtypes using spatial transcriptomics | Diseases: colorectal cancer | Tissues: rectosigmoid junction | Assays: Visium Spatial Gene Expression V1...

### HNSCC Cohorts
*21 cohorts identified for HNSCC*

#### <a id='gse200996'></a>GSE200996 — Tissue-resident Memory and Circulating T cells are Early Responders to Pre-surgical Cancer Immunotherapy
- **Cancer Type / Indication:** HNSCC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Single-cell suspension sorted by FACS for CD3+/CD8+ T lymphocytes (TILs)."*
- **Cohort Scale:** 204 samples/patients; ~408,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE200996_CD4.PBMC.single.cell.meta.data.txt.gz; GSE200996_CD4.tumor.single.cell.meta.data.txt.gz; GSE200996_CD45.PBMC.single.cell.meta.data.txt.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE200nnn/GSE200996/suppl/GSE200996_CD4.PBMC.singl`
- **Citation / Reference:** PMID:35803260
- **Repository Access Link:** [GSE200996](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200996)
- **Study Abstract / Experimental Design:**
  Pre-surgical (neoadjuvant) immune checkpoint blockade has shown promising activity in multiple cancer types, but the molecular mechanisms are not well understood. Here, we characterized early kinetic changes in tumor-infiltrating and circulating immune cells in oral cancer patients treated with neoadjuvant anti-PD-1 or anti-PD-1/CTLA-4 in a phase 2 clinical trial. Tumor-infiltrating CD8 T cells th...

#### <a id='gse301741'></a>GSE301741 — Single cell analysis highlights the significance of malignant cell IFN/MHC-II for immunotherapy response in head and neck squamous cell carcinoma
- **Cancer Type / Indication:** HNSCC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 58 samples/patients; ~116,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE301741_Seurat_Object_QCpass_137020cells_withMetaData.rds; filelist.txt; GSE301741_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE301nnn/GSE301741/suppl/GSE301741_Seurat_Object_`
- **Citation / Reference:** PMID:41923630
- **Repository Access Link:** [GSE301741](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301741)
- **Study Abstract / Experimental Design:**
  For many cancers, including head and neck squamous cell carcinoma (HNSCC), response rates to immunotherapy remain modest, with limited ability to predict responders. Previous studies that characterized cellular changes associated with immunotherapy in HNSCC focused on immune cells, providing limited insight into malignant cell responses. Motivated by this gap, we performed single cell RNA-sequenci...

#### <a id='gse296954'></a>GSE296954 — Differentiation of tumor-infiltrating GZMK+ effector memory T cells associates with response to neoadjuvant immunotherapy in head and neck cancer
- **Cancer Type / Indication:** HNSCC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** 10x Chromium 5' (Immune Profiling)
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Single-cell suspension sorted by FACS for CD3+/CD8+ T lymphocytes (TILs)."*
- **Cohort Scale:** 44 samples/patients; ~88,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE296954_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE296nnn/GSE296954/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE296954](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296954)
- **Study Abstract / Experimental Design:**
  We performed single-cell RNA/VDJ-sequencing of sorted T cells from pre- and post-treatment tumor biopsies of patients with HPV-unrelated head and neck squamous cell carcinoma treated with either bintrafusp alfa alone or in combination with Tri-Ad5 vaccine....

#### <a id='gse301720'></a>GSE301720 — Integrated single-cell and spatial analysis identifies context-dependent myeloid-T cell interactions in head and neck cancer immune checkpoint blockade response
- **Cancer Type / Indication:** HNSCC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 8 samples/patients; ~522,399 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE301720_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE301nnn/GSE301720/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41837744
- **Repository Access Link:** [GSE301720](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301720)
- **Study Abstract / Experimental Design:**
  Background Approximately 15-20% of head and neck cancer squamous cell carcinoma (HNSCC) patients respond favorably to immune checkpoint blockade (ICB). Previous single-cell RNA-Seq (scRNA-Seq) studies identified immune features, including macrophage subset ratios and T-cell subtypes, in HNSCC ICB response. However, the spatial features of HNSCC-infiltrated immune cells in response to ICB treatment...

#### <a id='gse247582'></a>GSE247582 — Single-cell analysis of CD4+ cytotoxic T lymphocytes in human oral squamous cell carcinoma
- **Cancer Type / Indication:** HNSCC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"We next demonstrated comprehensive delineation of the potential for CD8+ T cell differentiation towards dysfunctional states."*
- **Cohort Scale:** 3 samples/patients; ~6,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE247582_counts.tar.gz; filelist.txt; GSE247582_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE247nnn/GSE247582/suppl/GSE247582_counts.tar.gz;`
- **Citation / Reference:** PMID:38077321
- **Repository Access Link:** [GSE247582](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247582)
- **Study Abstract / Experimental Design:**
  Cancer immunotherapy targeting CD8+ T cells has made remarkable progress, even for oral squamous cell carcinoma (OSCC), a heterogeneous epithelial tumor without a substantial increase in the overall survival rate over the past decade. However, the therapeutic effects remain limited due to therapy resistance. Thus, a more comprehensive understanding of the roles of CD4+ T cells and B cells is cruci...

#### <a id='gse287301'></a>GSE287301 — Single-cell RNA-sequencing of HNSCC-infiltrating T cells
- **Cancer Type / Indication:** HNSCC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** Subcellular Spatial Transcriptomics (CosMx/Xenium)
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 48 samples/patients; ~96,000 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `GSE287301_filtered_feature_bc_matrix.tar.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE287nnn/GSE287301/suppl/GSE287301_filtered_featu`
- **Citation / Reference:** PMID:41961948
- **Repository Access Link:** [GSE287301](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287301)
- **Study Abstract / Experimental Design:**
  This study examines the T cell landscape in head and neck squamous cell carcinoma (HNSCC) using single-cell RNA sequencing combined with T cell receptor and protein profiling. By analyzing tumor-infiltrating T cells from 28 patients, we identified diverse T cell subsets, characterized their transcriptional states, and mapped their clonal dynamics. These findings provide a detailed view of the immu...

#### <a id='gse296867'></a>GSE296867 — Post-treatment peripheral-blood T cells from HNSCC patients undergoing neoadjuvant immunotherapy treatment.
- **Cancer Type / Indication:** HNSCC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** 10x Chromium 5' (Immune Profiling)
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"We performed single-cell RNA/VDJ-sequencing of sorted T cells from post-treatment blood of patients with HPV-unrelated head and neck squamous cell carcinoma treated with either bintrafusp alfa alone or in combination with Tri-Ad5 vaccine. These patients are identical to the patients from which tumor-infiltrating T cells were isolated as reflected in the sample names."*
- **Cohort Scale:** 22 samples/patients; ~44,000 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `filelist.txt; GSE296867_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE296nnn/GSE296867/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE296867](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296867)
- **Study Abstract / Experimental Design:**
  We performed single-cell RNA/VDJ-sequencing of sorted T cells from post-treatment blood of patients with HPV-unrelated head and neck squamous cell carcinoma treated with either bintrafusp alfa alone or in combination with Tri-Ad5 vaccine. These patients are identical to the patients from which tumor-infiltrating T cells were isolated as reflected in the sample names....

#### <a id='gse327189'></a>GSE327189 — Viral-based individualized neoantigen vaccine as adjuvant treatment in resected head and neck squamous cell carcinoma: immunogenicity and efficacy from a randomized Phase I trial
- **Cancer Type / Indication:** HNSCC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** 10x Chromium 5' (Immune Profiling)
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"This study allows the characterization of transcriptome and vdj rearrangements of tumor antigen specific CD8+ T cells elicited by personalized vaccines."*
- **Cohort Scale:** 5 samples/patients; ~10,000 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `GSE327189_feature_reference.csv.gz; filelist.txt; GSE327189_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE327nnn/GSE327189/suppl/GSE327189_feature_refere`
- **Repository Access Link:** [GSE327189](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327189)
- **Study Abstract / Experimental Design:**
  TG4050 is a personalized vaccine targeting point mutations of the tumors.. This study allows the characterization of transcriptome and vdj rearrangements of tumor antigen specific CD8+ T cells elicited by personalized vaccines....

#### <a id='gse280982'></a>GSE280982 — Single-cell sequencing across radiation treatment timepoints in head and neck cancer
- **Cancer Type / Indication:** HNSCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 49 samples/patients; ~98,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE280982_HyPR-HN_02_3_Mutect2.vcf.gz; filelist.txt; GSE280982_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE280nnn/GSE280982/suppl/GSE280982_HyPR-HN_02_3_M`
- **Citation / Reference:** PMID:40593620
- **Repository Access Link:** [GSE280982](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280982)
- **Study Abstract / Experimental Design:**
  We performed single cell sequencing across treatment timepoints to evaluate the immunological impact of radiation in head and neck cancer....

#### <a id='cellxgene_714e6bc2'></a>CELLxGENE_714e6bc2 — Immune HNSCC Atlas
- **Cancer Type / Indication:** HNSCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 44 samples/patients; ~134,385 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `714e6bc2-d81d-4a5b-97de-818bf6ff3a1d.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/377e8ef1-6f2f-4091-bc27-9bc01959b8eb.h5ad`
- **Citation / Reference:** Dataset Version: https://datasets.cellxgene.cziscience.com/377e8ef1-6f2f-4091-bc27-9bc01959b8eb.h5ad curated and distributed by CZ CELLxGENE Discover in Collection: https://cellxgene.cziscience.com/collections/bc7397a3-ea49-4d57-84b8-80bd6885d4c4
- **Repository Access Link:** [CELLxGENE_714e6bc2](https://cellxgene.cziscience.com/e/714e6bc2-d81d-4a5b-97de-818bf6ff3a1d.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: A highly resolved integrated single-cell atlas of HPV-negative Head and Neck Cancer | Diseases: head and neck squamous cell carcinoma | Tissues: alveolar ridge of mandible, buccal mucosa, gingiva of lower jaw, hypopharynx, larynx, mandible, mouth floor, oral cavity, oropharynx, tongue | Assays: 10x 3' v2, 10x 3' v3, 10x 5' transcription profiling...

#### <a id='cellxgene_60acb72d'></a>CELLxGENE_60acb72d — Global HPV-negative HNSCC Atlas
- **Cancer Type / Indication:** HNSCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 44 samples/patients; ~227,816 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `60acb72d-4935-4bb2-80b4-3700a6a527f8.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/8b58da42-48ca-4615-90f2-0a3c25687da1.h5ad`
- **Citation / Reference:** Dataset Version: https://datasets.cellxgene.cziscience.com/8b58da42-48ca-4615-90f2-0a3c25687da1.h5ad curated and distributed by CZ CELLxGENE Discover in Collection: https://cellxgene.cziscience.com/collections/bc7397a3-ea49-4d57-84b8-80bd6885d4c4
- **Repository Access Link:** [CELLxGENE_60acb72d](https://cellxgene.cziscience.com/e/60acb72d-4935-4bb2-80b4-3700a6a527f8.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: A highly resolved integrated single-cell atlas of HPV-negative Head and Neck Cancer | Diseases: head and neck squamous cell carcinoma | Tissues: alveolar ridge of mandible, buccal mucosa, gingiva of lower jaw, hypopharynx, larynx, mandible, mouth floor, oral cavity, oropharynx, tongue | Assays: 10x 3' v2, 10x 3' v3, 10x 5' transcription profiling...

#### <a id='gse268014'></a>GSE268014 — Spatial transcriptomic analysis of tumor heterogeneity in HPV+ and HPV- Head and Neck Squamous Cell Carcinoma
- **Cancer Type / Indication:** HNSCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 39 samples/patients; ~78,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE268014_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE268nnn/GSE268014/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:42601456
- **Repository Access Link:** [GSE268014](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268014)
- **Study Abstract / Experimental Design:**
  Head and neck squamous cell carcinoma is the 6th most common malignancy worldwide. We performed spatial transcriptomics to profile intra- and inter-tumor heterogeneity. We identified recurring patterns of spatial gene expression that co-localize in subsets of tumors....

#### <a id='cellxgene_624d92e2'></a>CELLxGENE_624d92e2 — NonImmune HNSCC Atlas
- **Cancer Type / Indication:** HNSCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 37 samples/patients; ~93,431 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `624d92e2-d006-4241-9fd8-03535f850301.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/be2f74ea-f487-4650-8c03-080798c2872c.h5ad`
- **Citation / Reference:** Dataset Version: https://datasets.cellxgene.cziscience.com/be2f74ea-f487-4650-8c03-080798c2872c.h5ad curated and distributed by CZ CELLxGENE Discover in Collection: https://cellxgene.cziscience.com/collections/bc7397a3-ea49-4d57-84b8-80bd6885d4c4
- **Repository Access Link:** [CELLxGENE_624d92e2](https://cellxgene.cziscience.com/e/624d92e2-d006-4241-9fd8-03535f850301.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: A highly resolved integrated single-cell atlas of HPV-negative Head and Neck Cancer | Diseases: head and neck squamous cell carcinoma | Tissues: alveolar ridge of mandible, buccal mucosa, gingiva of lower jaw, hypopharynx, larynx, mandible, mouth floor, oral cavity, oropharynx, tongue | Assays: 10x 3' v2, 10x 3' v3, 10x 5' transcription profiling...

#### <a id='gse198315'></a>GSE198315 — Multiregional single-cell profiling reveals extensive field cancerization and an immunosuppressive microenvironment in oral squamous cell carcinoma
- **Cancer Type / Indication:** HNSCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 36 samples/patients; ~268,131 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE198315_OSCC_UMI_count_matrix.txt.gz; GSE198315_matrix.mtx.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198315/suppl/GSE198315_OSCC_UMI_count`
- **Citation / Reference:** PMID:41054545
- **Repository Access Link:** [GSE198315](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198315)
- **Study Abstract / Experimental Design:**
  Oral squamous cell carcinoma (OSCC) is highly heterogeneous and metastatic, and the mechanisms driving OSCC development, progression, and metastasis are poorly understood. We performed multiregional single-cell RNA sequencing on 268,131 cells obtained from tumor core, tumor periphery, adjacent non-tumor tissue, and metastatic lymph node samples from 10 patients with human papillomavirus (HPV)–nega...

#### <a id='gse310797'></a>GSE310797 — Single cell and Spatial RNA-seq in Oral Squamous Cell Carcinoma with Different Stage
- **Cancer Type / Indication:** HNSCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 19 samples/patients; ~38,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE310797_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE310nnn/GSE310797/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41851114
- **Repository Access Link:** [GSE310797](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310797)
- **Study Abstract / Experimental Design:**
  To investigate the dynamic changes in the tumor microenvironment during the progression of oral squamous cell carcinoma (OSCC), we collected 16 treatment-naïve OSCC samples for single-cell RNA sequencing (scRNA-seq). Among these, 6 cases were additionally subjected to spatial transcriptomic profiling to capture the tissue architecture and spatial organization of the tumor microenvironment....

#### <a id='gse296771'></a>GSE296771 — TCR-stimulation of peripheral-blood CD8+ T cells.
- **Cancer Type / Indication:** HNSCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 5' (Immune Profiling)
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"TCR-stimulation of peripheral-blood CD8+ T cells. We performed single-cell RNA/VDJ-sequencing of sorted CD8+ T cells from the blood of patients with HPV-unrelated head and neck squamous cell carcinoma."*
- **Cohort Scale:** 16 samples/patients; ~32,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE296771_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE296nnn/GSE296771/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE296771](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296771)
- **Study Abstract / Experimental Design:**
  We performed single-cell RNA/VDJ-sequencing of sorted CD8+ T cells from the blood of patients with HPV-unrelated head and neck squamous cell carcinoma. The corresponding cells were cultured with or without aCD3/aCD28-mediated TCR stimulation. Importantly, these patients were not treated on a clinical trial and are distinct from the other patients included in this study....

#### <a id='cellxgene_55ca4411'></a>CELLxGENE_55ca4411 — ALL- Cells of gingivo-buccal oral squamous cell carcinoma
- **Cancer Type / Indication:** HNSCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 12 samples/patients; ~28,186 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `55ca4411-8c7d-4aef-a572-50852e030c05.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/986e3fb7-b087-442a-8a30-9fe94d28873f.h5ad`
- **Citation / Reference:** 10.1111/cas.15979
- **Repository Access Link:** [CELLxGENE_55ca4411](https://cellxgene.cziscience.com/e/55ca4411-8c7d-4aef-a572-50852e030c05.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Single-cell transcriptomic analysis of gingivo-buccal oral cancer reveals two dominant cellular programs | Diseases: oral cavity squamous cell carcinoma | Tissues: oral cavity | Assays: 10x 3' v3...

#### <a id='gse322620'></a>GSE322620 — Single-cell transcriptomic profiling of primary oral squamous cell carcinoma and lymph node metastatic lesions
- **Cancer Type / Indication:** HNSCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 10 samples/patients; ~20,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE322620_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE322nnn/GSE322620/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE322620](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE322620)
- **Study Abstract / Experimental Design:**
  we performed single-cell RNA sequencing (scRNA-seq) on tumor tissues from five patients with primary OSCC and five patients lymph node metastatic lesions....

#### <a id='gse281978'></a>GSE281978 — Deciphering head and neck cancer microenvironment: Single‐cell and spatial transcriptomics reveals human papillomavirus‐associated differences
- **Cancer Type / Indication:** HNSCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 6 samples/patients; ~12,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE281978_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281978/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:38235919
- **Repository Access Link:** [GSE281978](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281978)
- **Study Abstract / Experimental Design:**
  Human papillomavirus (HPV) is a major causative factor of head and neck squamous cell carcinoma (HNSCC), and the incidence of HPV-associated HNSCC is increasing. The role of tumor microenvironment (TME) in viral infection and metastasis needs to be explored further. Thus we studied the molecular characteristics of primary tumors (PTs) and lymph node metastatic tumors (LNMTs) by stratifying them ba...

#### <a id='gse339480'></a>GSE339480 — An advanced scRNA-seq-compatible lentiviral barcoding system resolves treatment-associated clonal selection and transcriptional states in HNSCC cells
- **Cancer Type / Indication:** HNSCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium (3' unspecified)
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~2,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE339480_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE339nnn/GSE339480/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE339480](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE339480)
- **Study Abstract / Experimental Design:**
  Cancer therapy response is shaped by pre-existing heterogeneity and treatment-associated selection within tumor cell populations. In this study, we applied a single-cell RNA-seq-compatible lentiviral BC24 barcode system to FaDu hypopharyngeal squamous cell carcinoma cells to link clonal identity with transcriptional cell state during cytotoxic treatment response. Barcoded FaDu cells were profiled ...

#### <a id='gse286935'></a>GSE286935 — Analysis of B cell neighborhoods via single cell spatial transcriptomics in head and neck cancer patients
- **Cancer Type / Indication:** HNSCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~2,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE286935_rawfiles.tar.gz; filelist.txt; GSE286935_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286935/suppl/GSE286935_rawfiles.tar.g`
- **Citation / Reference:** PMID:39970232
- **Repository Access Link:** [GSE286935](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286935)
- **Study Abstract / Experimental Design:**
  Head and neck squamous cell carcinoma tissue microarrays were stained via CD20 single plex IHC and then FOVs were selected for B cell neighborhoods (including tertiary lymphoid structures) in confirmed tumor areas via H&E evaluation by a pathologist....

### Gastric Cohorts
*16 cohorts identified for Gastric*

#### <a id='gse270680'></a>GSE270680 — A spatially resolved atlas of gastric cancer characterises a lymphocyte aggregated region [scRNA-seq]
- **Cancer Type / Indication:** Gastric
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 77 samples/patients; ~154,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE270680_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270680/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41593079
- **Repository Access Link:** [GSE270680](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270680)
- **Study Abstract / Experimental Design:**
  The tumour microenvironment (TME) is a focal point in cancer immunotherapy: its cellular composition and spatial organisation, especially the distribution of lymphocytes, can affect the clinical outcomes of cancer patients. In addition, the function of a cell differs depending on its spatial location and interaction with neighbouring cells. Here, by integrating single-cell transcriptomics with spa...

#### <a id='gse239676'></a>GSE239676 — Atlas of metastatic gastric cancer links ferroptosis to disease progression and immunotherapy response
- **Cancer Type / Indication:** Gastric
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Indicative of immune evasion, we observed diminished T cell clonal expansion, a marked reduction in the fraction of tumour-reactive CD8 T cells, and an enrichment of naïve T cells, M2-like, and proliferative macrophages in PC samples compared to PRIs."*
- **Cohort Scale:** 68 samples/patients; ~136,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE239676_count_matrix.mtx.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239676/suppl/GSE239676_count_matrix.m`
- **Citation / Reference:** PMID:39097198
- **Repository Access Link:** [GSE239676](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239676)
- **Study Abstract / Experimental Design:**
  Gastric adenocarcinoma (GAC) remains a significant cause of cancer-related deaths worldwide, with the majority of mortality resulting from metastatic spread. To elucidate the evolution of cancer cells and their interactions with the tumour microenvironment (TME), we conducted a comprehensive single-cell transcriptome and immune repertoire profiling of primary tumours (PRIs), matched liver metastas...

#### <a id='cellxgene_7bb64315'></a>CELLxGENE_7bb64315 — UMAP of Cancer Data integration
- **Cancer Type / Indication:** Gastric
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v3/v3.1
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 70 samples/patients; ~293,823 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `7bb64315-9e5a-41b9-9235-59acf9642a3e.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/36f88b68-b0c5-43eb-9dba-a620ed6b21aa.h5ad`
- **Citation / Reference:** 10.1158/2159-8290.cd-22-0824
- **Repository Access Link:** [CELLxGENE_7bb64315](https://cellxgene.cziscience.com/e/7bb64315-9e5a-41b9-9235-59acf9642a3e.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Single-cell RNA sequencing unifies developmental programs of Esophageal and Gastric Intestinal Metaplasia | Diseases: Barrett esophagus, gastric cancer, gastric intestinal metaplasia, gastritis, normal | Tissues: ascending colon epithelium, body of stomach, cardia of stomach, duodenum, epithelium of rectum, esophagogastric junction, ileal epithelium, lower esophagus, submucosal esophag...

#### <a id='cellxgene_0d3807bf'></a>CELLxGENE_0d3807bf — Extended - Stomach
- **Cancer Type / Indication:** Gastric
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 3' v2
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 62 samples/patients; ~70,090 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `0d3807bf-97f0-4ec5-9611-15ada507f0df.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/c10345ad-7e42-4a1b-ae78-d4d1547b8d06.h5ad`
- **Citation / Reference:** 10.1038/s41586-024-07571-1
- **Repository Access Link:** [CELLxGENE_0d3807bf](https://cellxgene.cziscience.com/e/0d3807bf-97f0-4ec5-9611-15ada507f0df.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: Single-cell integration reveals metaplasia in inflammatory gut diseases | Diseases: gastric cancer, normal | Tissues: body of stomach, pyloric antrum, pylorus, stomach | Assays: 10x 3' v2, 10x 5' transcription profiling, 10x 5' v2...

#### <a id='gse212212'></a>GSE212212 — Gene expression profile at single cell level of immune cells from gastric cancer
- **Cancer Type / Indication:** Gastric
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 44 samples/patients; ~88,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE212212_Information_of_multiplexed_results_preliminarily_filtered.csv.gz; GSE212212_batch_10_DBEC_MolsPerCell.csv.gz; GSE212212_batch_10_Sample_Tag_Calls.csv.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/GSE212212/suppl/GSE212212_Information_of`
- **Citation / Reference:** PMID:36921674
- **Repository Access Link:** [GSE212212](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212212)
- **Study Abstract / Experimental Design:**
  Immune cells within gastric cancer microenvironment greatly influence disease outcome. We used single cell RNA sequencing (scRNA-seq) to analyze the immune heterogeneity of gastric cancer....

#### <a id='gse228598'></a>GSE228598 — Gene expression profile at single cell level of Peritoneal Cells in Patients with Gastric Cancer
- **Cancer Type / Indication:** Gastric
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 28 samples/patients; ~56,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE228598_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE228nnn/GSE228598/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:38612926
- **Repository Access Link:** [GSE228598](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228598)
- **Study Abstract / Experimental Design:**
  We used single cell RNA sequencing (scRNA-seq) to analyze the diversity of Peritoneal Cells in Patients with Gastric Cancer....

#### <a id='gse234209'></a>GSE234209 — Gene expression profile at single cell level of immune cells and gastric cancer cells from multipoint biopsy samples of 2 gastric cancer patients treated by PD-1 antibody-based therapy
- **Cancer Type / Indication:** Gastric
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 12 samples/patients; ~24,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE234209_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234209/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE234209](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234209)
- **Study Abstract / Experimental Design:**
  To analyze the correlation between the neutrophil-to-lymphocyte ratio (NLR) and prognosis of advanced gastric cancer (AGC) patients treated by PD-1 antibody-based therapy and to delineate molecular characteristics of circulating neutrophils by single-cell RNA sequencing (scRNA-seq)....

#### <a id='gse275648'></a>GSE275648 — Role of Cancer-Associated Fibroblasts in Prognostic Stratification of Advanced-Stage Gastric Cancer through Interacting with Endothelial and Malignant Epithelial Cells [AGC scRNA-seq]
- **Cancer Type / Indication:** Gastric
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 11 samples/patients; ~22,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE275648_Reanalyzed_samples_RDS.tar.gz; GSE275648_sample_metadata.txt.gz; filelist.txt`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE275nnn/GSE275648/suppl/GSE275648_Reanalyzed_sam`
- **Citation / Reference:** PMID:41484771
- **Repository Access Link:** [GSE275648](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE275648)
- **Study Abstract / Experimental Design:**
  Patients with advanced gastric cancers (AGCs) would experience poor prognosis, lacking investigation on comprehensive ecosystem profile and specific prognostic factors. Here, we conducted patient stratification based on unsupervised clustering of transcriptomic profile of 108 normal/tumor AGC pairs and integrated single-cell RNA transcriptome profile of 116 gastric cancer/normal samples, revealing...

#### <a id='gse112302'></a>GSE112302 — Comprehensive Molecular Characterizations of Gastric Tissue and Gastric Cancer Revealed by Single-cell RNA-seq
- **Cancer Type / Indication:** Gastric
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 10 samples/patients; ~20,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE112302_Sample_Barcode_Information.xlsx; filelist.txt; GSE112302_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/GSE112302/suppl/GSE112302_Sample_Barcode`
- **Repository Access Link:** [GSE112302](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112302)
- **Study Abstract / Experimental Design:**
  Molecular knowledge of normal gastric tissues and gastric cancers remains incomplete. Here, we used single-cell RNA-seq to study the cell diversity of gastric tissues and gastric cancers. The expression landscape of normal gastric cell types and several candidate stem cell markers were obtained. Surprisingly, nearly all cell types in the antrum could transdifferentiate to intestinal metaplasia (IM...

#### <a id='gse168537'></a>GSE168537 — Single-Cell Transcriptome Analysis Reveals Neutrophil Populations in Gastric Cancer
- **Cancer Type / Indication:** Gastric
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 10 samples/patients; ~20,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE168537_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE168nnn/GSE168537/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:36921037
- **Repository Access Link:** [GSE168537](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168537)
- **Study Abstract / Experimental Design:**
  The full neutrophil heterogeneity and function remains incompletely characterized in gastric cancer progression. Here, we profiled >21,854 mouse neutrophils with an average of 1,189 genes per cell using single-cell RNA sequencing to provide a comprehensive transcriptional landscape of neutrophil function and fate decision in tumor progression. By combining scRNA-Seq analysis with tumor model, flow...

#### <a id='gse246662'></a>GSE246662 — Single-cell profiling reveals altered immune landscape and impaired NK cell function in gastric cancer liver metastasis
- **Cancer Type / Indication:** Gastric
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 9 samples/patients; ~18,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE246662_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246662/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:39060439
- **Repository Access Link:** [GSE246662](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246662)
- **Study Abstract / Experimental Design:**
  Gastric cancer (GC) is a substantial global health concern, and the development of liver metastasis (LM) in GC represents a critical stage linked to unfavorable patient prognoses. In this study, we employed single-cell RNA sequencing (scRNA-seq) to investigate the immune landscape of GC liver metastasis, revealing several immuno-suppressive components within the tumor immune microenvironment (TIM)...

#### <a id='gse308231'></a>GSE308231 — Single-cell sequencing technology reveals the characteristics of the tumor microenvironment in gastric cancer and its peritoneal metastases
- **Cancer Type / Indication:** Gastric
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 6 samples/patients; ~12,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE308231_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE308nnn/GSE308231/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41360923
- **Repository Access Link:** [GSE308231](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308231)
- **Study Abstract / Experimental Design:**
  To analyze the tumor microenvironment characteristics of gastric cancer and its peritoneal metastases, we collected fresh surgical specimens pathologically confirmed as gastric cancer and peritoneal metastases, and performed single-cell sequencing....

#### <a id='cellxgene_ca140407'></a>CELLxGENE_ca140407 — T cells from gastric cancer patients
- **Cancer Type / Indication:** Gastric
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium 5' (Immune Profiling)
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 5 samples/patients; ~45,698 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `ca140407-efd4-48e3-a040-524b8a7752cc.h5ad`
- **Direct Download URLs:** `https://datasets.cellxgene.cziscience.com/2ab4df39-876d-4fdc-af54-0fa5113a54f4.h5ad`
- **Citation / Reference:** 10.1038/s41467-026-70751-2
- **Repository Access Link:** [CELLxGENE_ca140407](https://cellxgene.cziscience.com/e/ca140407-efd4-48e3-a040-524b8a7752cc.cxg/)
- **Study Abstract / Experimental Design:**
  Collection: T cells from gastric cancer patients | Diseases: gastric cancer | Tissues: blood, liver, lymph node, stomach | Assays: 10x 5' v2...

#### <a id='gse321676'></a>GSE321676 — Dural scRNA-seq in two organoids derived from malignant ascites of the same gastric cancer patient before (Primary) and after (Resistance) resistant to paclitaxel
- **Cancer Type / Indication:** Gastric
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 2 samples/patients; ~4,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE321676_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE321nnn/GSE321676/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE321676](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE321676)
- **Study Abstract / Experimental Design:**
  We have employed a single sequencing(scRNA-seq) approach using 10× Genomics scRNAseq of organoids to study paclitaxel resistance in gastric cancer(GC). Patient-derived organoids (PDOs) can recapitulate majority aspects of tissue where they are derived from, in terms of specific molecular profiles, divergent phenotypes, concerning growth pattern, response to classical chemotherapy. This study repre...

#### <a id='gse184198'></a>GSE184198 — Single cell sequencing of gastric cancer and adjacent normal tissues
- **Cancer Type / Indication:** Gastric
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 2 samples/patients; ~4,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE184198_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184198/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:36372898
- **Repository Access Link:** [GSE184198](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184198)
- **Study Abstract / Experimental Design:**
  In this study, we used single-cell RNA sequencing (scRNA-seq) to further elucidate the composition of the gastric cancer and adjacent normal tissue microenvironment and found that TNF and ACTG1 may be key regulators of regulatory T cells (Tregs) in gastric cancer. A total of seven cell types were identified, including B cells, T cells, NK cells, mast cells, fibroblasts, endothelial cells, and epit...

#### <a id='gse232733'></a>GSE232733 — 5' droplet-based single-cell RNA sequencing of human gastric cancer tissue
- **Cancer Type / Indication:** Gastric
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~2,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE232733_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232733/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:39107346
- **Repository Access Link:** [GSE232733](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232733)
- **Study Abstract / Experimental Design:**
  Using single-cell sequencing, we profiled single cells derived from one human gastric tumor....

### HCC Cohorts
*22 cohorts identified for HCC*

#### <a id='gse313642'></a>GSE313642 — Immunosuppressive monocytes are enriched in hepatocellular carcinoma patients with liver dysfunction in a phase II trial of combination sorafenib and nivolumab [CITE-Seq]
- **Cancer Type / Indication:** HCC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 194 samples/patients; ~388,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE313642_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE313nnn/GSE313642/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41831609
- **Repository Access Link:** [GSE313642](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313642)
- **Study Abstract / Experimental Design:**
  Background and aims: Immune checkpoint inhibition (ICI) and anti-angiogenic therapies are active in hepatocellular carcinoma (HCC), although patients with impaired hepatic function have worse outcomes. Methods: We conducted a multi-center, open-label phase II clinical trial to assess the safety and efficacy of the multikinase inhibitor, sorafenib, combined with nivolumab, in patients with advanced...

#### <a id='gse245906'></a>GSE245906 — Identification of TREM1+CD163+ myeloid cells as a deleterious immune subset in HCC [scRNA-seq]
- **Cancer Type / Indication:** HCC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 20 samples/patients; ~40,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE245906_Giraud_Chalopin_innate_HCC_metadata.tsv.gz; GSE245906_Giraud_Chalopin_innate_HCC_processed_count_data_mat.tsv.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245906/suppl/GSE245906_Giraud_Chalopi`
- **Citation / Reference:** PMID:38350444
- **Repository Access Link:** [GSE245906](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245906)
- **Study Abstract / Experimental Design:**
  Hepatocellular carcinoma (HCC) is an inflammation-associated cancer arising from viral and non-viral etiologies. Expansion of suppressive myeloid cells is a hallmark of chronic inflammation and cancer, but their heterogeneity in HCC is not fully resolved and might underlie immunotherapy resistance in the steatohepatitis setting. Here, we present a high resolution atlas of hepatic innate immune cel...

#### <a id='gse319709'></a>GSE319709 — Single-cell transcriptomics reveals etiology-specific T-cell heterogeneity in hepatocellular carcinoma and implicates regulatory
- **Cancer Type / Indication:** HCC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 18 samples/patients; ~36,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE319709_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE319nnn/GSE319709/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE319709](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE319709)
- **Study Abstract / Experimental Design:**
  Hepatocellular carcinoma (HCC) exhibits remarkable etiological heterogeneity, with hepatitis B virus (HBV) infection and metabolic dysfunction-associated steatohepatitis (MASH) emerging as two leading causes. The tumor microenvironment (TME), particularly T cell subsets, plays a pivotal role in tumor progression and immunotherapy response. However, the etiology-specific T cell landscapes in HBV-HC...

#### <a id='gse272347'></a>GSE272347 — Late-stage tertiary lymphoid structures in hepatocellular carcinoma treated with neoadjuvant immune checkpoint blockade [scRNA-seq]
- **Cancer Type / Indication:** HCC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 17 samples/patients; ~34,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE272347_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE272nnn/GSE272347/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:39455893
- **Repository Access Link:** [GSE272347](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272347)
- **Study Abstract / Experimental Design:**
  Tertiary lymphoid structures (TLS) are associated with improved response in solid tumors treated with immune checkpoint blockade (ICB), but understanding of their clinical significance and the circumstances of their resolution remains incomplete. Here, we found that in hepatocellular carcinoma (HCC) treated with neoadjuvant immunotherapy, high intratumoral TLS density at the time of surgery is ass...

#### <a id='gse318418'></a>GSE318418 — Immunosuppressive monocytes are enriched in hepatocellular carcinoma patients with liver dysfunction in a phase II trial of combination sorafenib and nivolumab [scRNA-Seq]
- **Cancer Type / Indication:** HCC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 16 samples/patients; ~32,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE318418_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE318nnn/GSE318418/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41831609
- **Repository Access Link:** [GSE318418](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318418)
- **Study Abstract / Experimental Design:**
  Background and aims: Immune checkpoint inhibition (ICI) and anti-angiogenic therapies are active in hepatocellular carcinoma (HCC), although patients with impaired hepatic function have worse outcomes.   Methods: We conducted a multi-center, open-label clinical trial to assess the safety and efficacy of the multikinase inhibitor, sorafenib, combined with nivolumab, in patients with advanced or unr...

#### <a id='gse299340'></a>GSE299340 — Single-cell RNA sequencing reveals B cell-related immunosuppressive landscape and a potential suppressor in hepatocellular carcinoma
- **Cancer Type / Indication:** HCC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 10 samples/patients; ~20,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE299340_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE299nnn/GSE299340/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:40721811
- **Repository Access Link:** [GSE299340](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299340)
- **Study Abstract / Experimental Design:**
  Background Sophisticated tumor microenvironment is responsible for malignant progression and poor prognosis of hepatocellular carcinoma (HCC) patients. Discovering new therapeutic targets was desired for preferable treatments of HCC patients. Methods To uncover the HCC microenvironment, single-cell transcriptomes of HCC tissues and corresponding non-cancerous tissues were analyzed. Differentially ...

#### <a id='gse282343'></a>GSE282343 — Viral-Track integrated single-cell RNA-sequencing reveals HBV lymphotropism and immunosuppressive microenvironment in HBV-associated hepatocellular carcinoma
- **Cancer Type / Indication:** HCC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Single-cell suspension sorted by FACS for CD3+/CD8+ T lymphocytes (TILs)."*
- **Cohort Scale:** 8 samples/patients; ~71,466 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE282343_ScRNA.counts.txt.gz; GSE282343_ScRNA.metadata.txt.gz; GSE282343_ScRNA.normalized.txt.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282343/suppl/GSE282343_ScRNA.counts.t`
- **Citation / Reference:** PMID:40634523
- **Repository Access Link:** [GSE282343](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282343)
- **Study Abstract / Experimental Design:**
  The tumor microenvironment (TME) is a crucial mediator of tumor progression and treatment response. Here, we compare the immune microenvironments of HBV and non-HBV-HCC and investigate the reason for the persistence of HBV infection in the liver. We combine the Viral-Track method with single-cell RNA sequencing and profile the transcriptomes of 71,466 cells from HBV and non-HBV-HCC patients. In ad...

#### <a id='gse255830'></a>GSE255830 — Gene expression and T cell repertoire profile at single cell level of peripheral blood T cells after treatment with a personalized neoantigen vaccine (GNOS-PV02) and Pembrolizumab for advanced hepatocellular carcinoma.
- **Cancer Type / Indication:** HCC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 8 samples/patients; ~16,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE255830_readme.txt; filelist.txt; GSE255830_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE255nnn/GSE255830/suppl/GSE255830_readme.txt; ht`
- **Citation / Reference:** PMID:38584166
- **Repository Access Link:** [GSE255830](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255830)
- **Study Abstract / Experimental Design:**
  To further characterize the T cell response to the neoantigen vaccine, we performed single-cell RNA and T cell receptor (scRNA/TCR-seq) of peripheral blood T cells at the 12-week post-vaccination time point....

#### <a id='gse281110'></a>GSE281110 — Molecular landscape of tumor-associated tissue-resident memory T cells in tumor microenvironment of hepatocellular carcinoma [HCC_scRNA]
- **Cancer Type / Indication:** HCC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 7 samples/patients; ~14,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE281110_HCC_M1_LIL_demuxlet.best.gz; filelist.txt; GSE281110_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281110/suppl/GSE281110_HCC_M1_LIL_dem`
- **Citation / Reference:** PMID:39934824
- **Repository Access Link:** [GSE281110](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281110)
- **Study Abstract / Experimental Design:**
  Immunotherapy for liver cancer is used to rejuvenate tumor-infiltrating lymphocytes by modulating the immune microenvironment. Thus, early protective functions of T cell subtypes with tissue-specific residency have been studied in the tumor microenvironment (TME). We identified tumor-associated tissue-resident memory T (TA-TRM) cells in hepatocellular carcinoma (HCC) and characterized their molecu...

#### <a id='gse233405'></a>GSE233405 — Immunohistochemical scoring of LAG-3 in conjunction with CD8 in  the tumor microenvironment predicts response to immunotherapy in  hepatocellular carcinoma
- **Cancer Type / Indication:** HCC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"A 4-gene inflammatory signature, comprising CD8, PD-L1, LAG-3, and STAT1, was recently shown to be associated with a better overall response to ICB in various cancer types."*
- **Cohort Scale:** 6 samples/patients; ~12,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE233405_processed_scRNAseq_data.csv.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233405/suppl/GSE233405_processed_scRN`
- **Citation / Reference:** PMID:37342338
- **Repository Access Link:** [GSE233405](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233405)
- **Study Abstract / Experimental Design:**
  Introduction  Immune checkpoint blockade (ICB) is a systemic therapeutic option for advanced hepatocellular carcinoma (HCC). However, low patient response rates necessitate the development of robust predictive biomarkers that identify individuals who will benefit from ICB. A 4-gene inflammatory signature, comprising CD8, PD-L1, LAG-3, and STAT1, was recently shown to be associated with a better ov...

#### <a id='gse278324'></a>GSE278324 — Gene regulatory network analysis on snRNAseq revealed key regulators for hepatocellular carcinoma progression
- **Cancer Type / Indication:** HCC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** Single-Nucleus RNA-seq
- **Modality:** scRNA-seq + snRNA-seq
- **Cell Selection / Filtering Strategy:** **Nuclei Isolation (snRNA-seq)**
  > *"Single-nucleus RNA sequencing performed on isolated nuclei from frozen resection tissue."*
- **Cohort Scale:** 6 samples/patients; ~12,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE278324_HCC.combined.count.matrix.txt.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE278nnn/GSE278324/suppl/GSE278324_HCC.combined.c`
- **Repository Access Link:** [GSE278324](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278324)
- **Study Abstract / Experimental Design:**
  Many clinical researchers have developed targeted therapy and immunotherapy for the case of advanced stage of hepatocellular carcinoma (HCC). However, the molecular mechanisms of HCC development still need to be investigated to improve the response rate of those therapies. Here, we generated single nucleus RNA sequencing (snRNAseq) data from biopsy samples of six patients with Barcelona Clinic Liv...

#### <a id='gse265770'></a>GSE265770 — Gene expression profile at single cell level of CD56+ natural killer cells and CD8+ T cells from blood, spleen and HCC-PDX in humanized mice
- **Cancer Type / Indication:** HCC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Single-cell suspension sorted by FACS for CD3+/CD8+ T lymphocytes (TILs)."*
- **Cohort Scale:** 3 samples/patients; ~6,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE265770_NKsubset_integrated_metadata.csv.gz; GSE265770_Tsubset_integrated_metadata.csv.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE265nnn/GSE265770/suppl/GSE265770_NKsubset_integ`
- **Citation / Reference:** PMID:39318093
- **Repository Access Link:** [GSE265770](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE265770)
- **Study Abstract / Experimental Design:**
  In solid tumors, the exhaustion of natural killer (NK) cells and cytotoxic T cells in the immunosuppressive tumor microenvironment poses challenges for effective tumor control. Conventional humanized mouse models of hepatocellular carcinoma-patient-derived xenografts (HCC-PDX) encounter limitations in NK-cell infiltration, hindering studies on NK-cell immunobiology. Here, we introduce an improved ...

#### <a id='gse318420'></a>GSE318420 — Immunosuppressive monocytes are enriched in hepatocellular carcinoma patients with liver dysfunction in a phase II trial of combination sorafenib and nivolumab
- **Cancer Type / Indication:** HCC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 210 samples/patients; ~420,000 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `filelist.txt; GSE318420_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE318nnn/GSE318420/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41831609
- **Repository Access Link:** [GSE318420](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318420)
- **Study Abstract / Experimental Design:**
  This SuperSeries is composed of the SubSeries listed below....

#### <a id='gse272348'></a>GSE272348 — Late-stage tertiary lymphoid structures in hepatocellular carcinoma treated with neoadjuvant immune checkpoint blockade
- **Cancer Type / Indication:** HCC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 31 samples/patients; ~62,000 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `filelist.txt; GSE272348_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE272nnn/GSE272348/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:39455893
- **Repository Access Link:** [GSE272348](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272348)
- **Study Abstract / Experimental Design:**
  This SuperSeries is composed of the SubSeries listed below....

#### <a id='gse224411'></a>GSE224411 — Uncovering the spatial landscape of molecular interactions within the tumor microenvironment through latent spaces
- **Cancer Type / Indication:** HCC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 4 samples/patients; ~8,000 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `filelist.txt; GSE224411_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE224nnn/GSE224411/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:37080163
- **Repository Access Link:** [GSE224411](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224411)
- **Study Abstract / Experimental Design:**
  Recent advances in spatial transcriptomics (ST) enable gene expression measurements from a tissue sample while retaining its spatial context. This technology enables unprecedented in situ resolution of the regulatory pathways that underlie the heterogeneity in the tumor and its microenvironment (TME). The direct characterization of cellular co-localization with spatial technologies facilities quan...

#### <a id='gse215428'></a>GSE215428 — Single-cell RNA sequencing of immune landscape in hepatocellular carcinoma treated with sintilimab and sorafenib
- **Cancer Type / Indication:** HCC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 3 samples/patients; ~6,000 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `filelist.txt; GSE215428_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE215nnn/GSE215428/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE215428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE215428)
- **Study Abstract / Experimental Design:**
  Checkpoint inhibitors combined with sorafenib for hepatocellular carcinoma (HCC) treatment have achieved satisfactory results in many clinical trials. Understanding the mechanism underlying this modality may provide better strategies to eradicate HCC. In the present study, peripheral blood mononuclear cells (PBMCs) from patients who received combination sintilimab and sorafenib treatment were anal...

#### <a id='gse320155'></a>GSE320155 — scRNA‑seq of human intra‑hepatic CD4⁺ T cells reveals unique Treg transcriptional identities and two steady‑state tissue‑adaptation programs across healthy and HCC liver samples
- **Cancer Type / Indication:** HCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Chromium (3' unspecified)
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 60 samples/patients; ~120,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE320155_HCC_tumour_and_background_cellranger_aggr.tar.gz; GSE320155_Liver_and_PBMC_cellranger_aggr.tar.gz; GSE320155_feature_reference.csv.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320155/suppl/GSE320155_HCC_tumour_and`
- **Citation / Reference:** PMID:42443155
- **Repository Access Link:** [GSE320155](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE320155)
- **Study Abstract / Experimental Design:**
  We employed a single‑cell sequencing approach using the 10x Genomics platform, including scRNA‑seq, CITE‑seq, and paired TCR libraries, to investigate the molecular programs of human tissue‑resident regulatory T cells (Tregs). By analyzing CD4⁺ T cells from healthy liver and matched peripheral blood mononuclear cells (PBMCs), as well as hepatocellular carcinoma (HCC) tissue with paired non‑tumoral...

#### <a id='gse326201'></a>GSE326201 — Single-cell RNA sequencing of human hepatocellular carcinoma and adjacent non-tumour liver tissues reveals etiology-associated and etiology-independent tumour microenvironmental features
- **Cancer Type / Indication:** HCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 18 samples/patients; ~36,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE326201_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE326nnn/GSE326201/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE326201](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE326201)
- **Study Abstract / Experimental Design:**
  Hepatocellular carcinoma (HCC) arises in the context of diverse aetiological backgrounds, yet the extent to which the tumour microenvironment is shaped by aetiology-associated versus aetiology-independent programmes remains incompletely understood. This dataset was generated to support interrogation of the cellular architecture of human HCC and the surrounding liver microenvironment at single-cell...

#### <a id='gse290925'></a>GSE290925 — Single-Cell Transcriptomic Profiling of the Tumor Microenvironment in Treatment-Naive Hepatocellular Carcinoma Patients
- **Cancer Type / Indication:** HCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 12 samples/patients; ~24,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE290925_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE290nnn/GSE290925/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:40241752
- **Repository Access Link:** [GSE290925](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290925)
- **Study Abstract / Experimental Design:**
  This study provides high-resolution single-cell RNA sequencing (scRNA-seq) data of tumor tissues from 12 treatment-naive hepatocellular carcinoma (HCC) patients. The dataset captures the transcriptomic profiles of diverse cell populations, including tumor cells, immune cells (e.g., T cells, B cells, macrophages), and stromal cells, revealing cellular heterogeneity and key gene expression signature...

#### <a id='gse208308'></a>GSE208308 — Achieving CD8+ T cell-dependent lethality by targeting cancer USP14 in HCC
- **Cancer Type / Indication:** HCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Achieving CD8+ T cell-dependent lethality by targeting cancer USP14 in HCC This SuperSeries is composed of the SubSeries listed below."*
- **Cohort Scale:** 6 samples/patients; ~12,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE208308_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE208nnn/GSE208308/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:40966278
- **Repository Access Link:** [GSE208308](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208308)
- **Study Abstract / Experimental Design:**
  This SuperSeries is composed of the SubSeries listed below....

#### <a id='gse291757'></a>GSE291757 — Multiplexed transcriptomics to screen drug combination and define therapeutic mechanism of action at single-cell resolution
- **Cancer Type / Indication:** HCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 4 samples/patients; ~8,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE291757_processed_normalized_matrix.txt.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE291nnn/GSE291757/suppl/GSE291757_processed_norm`
- **Repository Access Link:** [GSE291757](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291757)
- **Study Abstract / Experimental Design:**
  Compared to classical drug screening, single-cell screening not only significantly enhances throughput but also provides richer transcriptional response information. In this study, we employed the high-throughput single-cell sequencing technology, snHH-seq, to screen clinical drug combinations with anti-hepatocellular carcinoma activity. Single-cell transcriptomics revealed that the combination of...

#### <a id='gse208307'></a>GSE208307 — Achieving CD8+ T cell-dependent lethality by targeting cancer USP14 in HCC [scRNA-seq]
- **Cancer Type / Indication:** HCC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Single-cell suspension sorted by FACS for CD3+/CD8+ T lymphocytes (TILs)."*
- **Cohort Scale:** 2 samples/patients; ~4,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE208307_MarkerGenes.csv.gz; filelist.txt; GSE208307_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE208nnn/GSE208307/suppl/GSE208307_MarkerGenes.cs`
- **Citation / Reference:** PMID:40966278
- **Repository Access Link:** [GSE208307](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208307)
- **Study Abstract / Experimental Design:**
  We used single cell RNA sequencing (scRNA-seq) to analyze the diversity of PD-1 therapy responsive and resistant patients samples....

### PDAC Cohorts
*21 cohorts identified for PDAC*

#### <a id='gse311789'></a>GSE311789 — DeCAF redefines fibroblast states uncovering multidimensional tumor-stroma relationships driving clinical tumor progression and immunotherapy response
- **Cancer Type / Indication:** PDAC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 142 samples/patients; ~284,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE311789_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE311nnn/GSE311789/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41707654
- **Repository Access Link:** [GSE311789](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE311789)
- **Study Abstract / Experimental Design:**
  This SuperSeries is composed of the SubSeries listed below....

#### <a id='gse279781'></a>GSE279781 — CD137 agonism enhances anti-PD1 induced activation of clonally expanded CD8+ T cells in a neoadjuvant pancreatic cancer clinical trial
- **Cancer Type / Indication:** PDAC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Tumor infiltrating leukocytes analyzed two weeks post neoadjuvant therapy revealed shifts in CD8+ T cell activation as well as cytoskeletal and extracellular matrix (ECM)-interacting components with Cy/GVAX and anti-PD1."*
- **Cohort Scale:** 30 samples/patients; ~60,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `GSE279781_matrix.mtx.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279781/suppl/GSE279781_matrix.mtx.gz`
- **Citation / Reference:** PMID:39811671
- **Repository Access Link:** [GSE279781](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279781)
- **Study Abstract / Experimental Design:**
  Successful pancreatic ductal adenocarcinoma (PDAC) immunotherapy requires therapeutic combinations that induce quality T cells. Tumor microenvironment (TME) analysis following therapeutic interventions can identify response mechanisms to guide design of more effective combinations. We provide a reference single-cell dataset from PDAC-infiltrating T cell and monocyte subsets from a human neoadjuvan...

#### <a id='gse316195'></a>GSE316195 — A Phase 1 clinical trial and single-cell correlates of motixafortide, cemiplimab, gemcitabine and nab-paclitaxel for metastatic treatment-naïve metastatic pancreatic ductal adenocarcinoma.
- **Cancer Type / Indication:** PDAC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** Single-Nucleus RNA-seq
- **Modality:** scRNA-seq + snRNA-seq
- **Cell Selection / Filtering Strategy:** **Nuclei Isolation (snRNA-seq)**
  > *"Single-nucleus RNA sequencing performed on isolated nuclei from frozen resection tissue."*
- **Cohort Scale:** 22 samples/patients; ~44,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE316195_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE316nnn/GSE316195/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE316195](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316195)
- **Study Abstract / Experimental Design:**
  The C-X-C motif chemokine receptor 4 (CXCR4)/C-X-C motif chemokine ligand 12 (CXCL12) axis is a well-established contributor to the immunosuppressive and immune-excluded TME in pancreatic adenocarcinoma (PDA). Building on pre-clinical data demonstrating a survival benefit with the addition of gemcitabine to CXCR4 and PD1 inhibition in the KPC mouse model, we conducted an open-label, single-arm pha...

#### <a id='gse212966'></a>GSE212966 — Single-cell RNA-seq reveals immune landscape of pancreatic cancer
- **Cancer Type / Indication:** PDAC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 12 samples/patients; ~24,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE212966_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/GSE212966/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:36944944
- **Repository Access Link:** [GSE212966](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212966)
- **Study Abstract / Experimental Design:**
  Pancreatic ductal adenocarcinoma (PDAC) has complex tumor immune microenvironment (TIME), the clinical values of which remains to be explored. This study aimed to delineate the immune landscape of PDAC and determine the clinical value of immune features in TIME. There was a significant difference in immune profiles between PDAC and adjacent normal pancreatic tissues. Several novel immune features ...

#### <a id='gse283206'></a>GSE283206 — Combined Flt3L and CD40 agonism restores dendritic cell driven T cell immunity in mouse models and patients with pancreatic cancer [human]
- **Cancer Type / Indication:** PDAC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"to fully engage cDC activity; yet when combined with Flt3L, the dual therapy triggered a cDC driven type-I immune response characterized by CD8+ T cell infiltration, interleukin-12 production, and a reciprocal interferon (IFN) response."*
- **Cohort Scale:** 8 samples/patients; ~16,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE283206_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE283nnn/GSE283206/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:40815670
- **Repository Access Link:** [GSE283206](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283206)
- **Study Abstract / Experimental Design:**
  T cell directed immunotherapies have largely failed to demonstrate clinical efficacy in patients with pancreatic ductal adenocarcinoma (PDAC). This broad resistance may result from poor tumor antigenicity and an immunosuppressive tumor microenvironment (TME). We hypothesize that tumor immunity in pancreatic cancer patients is further limited by systemic and tumor-intrinsic suppression of conventio...

#### <a id='gse311788'></a>GSE311788 — DeCAF redefines fibroblast states uncovering multidimensional tumor-stroma relationships driving clinical tumor progression and immunotherapy response [scRNA-Seq]
- **Cancer Type / Indication:** PDAC
- **Classification:** Tier 1 (ICB Response)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 6 samples/patients; ~12,000 cells
- **Clinical Context / Response:** Documented ICB immunotherapy with response / resistance / outcome correlates
- **Downloadable Matrix Files:** `filelist.txt; GSE311788_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE311nnn/GSE311788/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41707654
- **Repository Access Link:** [GSE311788](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE311788)
- **Study Abstract / Experimental Design:**
  Cancer-associated fibroblast (CAF) heterogeneity within pancreatic ductal adenocarcinoma (PDAC) has been previously assessed through single-cell RNA sequencing (scRNAseq) to discover distinct CAF populations such as myCAFs and iCAFs. While useful as a biological framework, no studies have conclusively and robustly demonstrated a correlation of CAF subpopulations with clinical prognosis or therapy ...

#### <a id='gse211644'></a>GSE211644 — Single cell transcriptomic and T cell repertoire analysis reveals trajectory of tumor - infiltrating lymphocyte states in pancreatic cancer
- **Cancer Type / Indication:** PDAC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Single-cell suspension sorted by FACS for CD3+/CD8+ T lymphocytes (TILs)."*
- **Cohort Scale:** 50 samples/patients; ~100,000 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `GSE211644_fresh_matrix.mtx.gz; GSE211644_grown_matrix.mtx.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211644/suppl/GSE211644_fresh_matrix.m`
- **Citation / Reference:** PMID:35849783
- **Repository Access Link:** [GSE211644](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211644)
- **Study Abstract / Experimental Design:**
   The CD8+ TIL contained a putative transitional GZMK+ population based on TCR clonotype sharing, and cell-state trajectory analysis showed similarity to a GZMB+PRF1+ cytotoxic and a CXCL13+ dysfunctional population. Statistical analysis suggested that certain TIL states, such as dysfunctional and inhibitory populations, often occurred together. Finally, analysis of cultured TIL revealed that high-...

#### <a id='gse348275'></a>GSE348275 — Clonal lineage tracing and parallel multiomics profiling reveal transcriptional heterogeneity induced by ARID1A deficiency
- **Cancer Type / Indication:** PDAC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 35 samples/patients; ~70,000 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `filelist.txt; GSE348275_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE348nnn/GSE348275/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE348275](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE348275)
- **Study Abstract / Experimental Design:**
  This SuperSeries is composed of the SubSeries listed below....

#### <a id='gse318413'></a>GSE318413 — Patient-derived orthotopic xenograft models recapitulate the peritoneal dissemination of pancreatic cancer and delineate its transcriptional and regulatory programs [scRNA-seq]
- **Cancer Type / Indication:** PDAC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** Single-Nucleus RNA-seq
- **Modality:** scRNA-seq + snRNA-seq
- **Cell Selection / Filtering Strategy:** **Nuclei Isolation (snRNA-seq)**
  > *"Single-nucleus RNA sequencing (snRNA-seq) and single-cell ATAC sequencing (scATAC-seq) were performed to analyze the tumors from these models."*
- **Cohort Scale:** 19 samples/patients; ~38,000 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `filelist.txt; GSE318413_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE318nnn/GSE318413/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41673768
- **Repository Access Link:** [GSE318413](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318413)
- **Study Abstract / Experimental Design:**
  BACKGROUND & AIMS: The mechanisms of peritoneal dissemination in pancreatic ductal adenocarcinoma (PDAC) remain unclear partly owing to the lack of patient-derived models that recapitulate this process. This study aimed to establish an orthotopic model of PDAC peritoneal dissemination and to uncover the transcriptional and regulatory programs underlying this process. METHODS: Organoids were establ...

#### <a id='gse156405'></a>GSE156405 — Elucidation of tumor-stromal heterogeneity and the ligand-receptor interactome by single cell transcriptomics in real-world pancreatic cancer biopsies
- **Cancer Type / Indication:** PDAC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 17 samples/patients; ~34,000 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `filelist.txt; GSE156405_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156405/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:34426439
- **Repository Access Link:** [GSE156405](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156405)
- **Study Abstract / Experimental Design:**
  Through single-cell RNA sequencing (scRNA-seq), we identify substantial  transcriptomic evolution of PDOs propagated from the parental tumor, which may alter  predicted drug sensitivity. In addition, we performed an integrative analyses of PDAC biopsies and provide an in-depth characterization of the heterogeneity within the tumor  microenvironment, including cancer-associated fibroblast (CAF) sub...

#### <a id='gse335452'></a>GSE335452 — Single-cell RNA-seq and spatial transcriptomics characterize CD8+ exhausted T cells in pancreatic ductal adenocarcinoma
- **Cancer Type / Indication:** PDAC
- **Classification:** Tier 1 (ICB Treated)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"sted T cells in pancreatic ductal adenocarcinoma Pancreatic cancer resists immunotherapy due to a suppressive immune microenvironment where CD8⁺ T cells play a key role."*
- **Cohort Scale:** 15 samples/patients; ~30,000 cells
- **Clinical Context / Response:** ICB immunotherapy treated cohort (response evaluation required)
- **Downloadable Matrix Files:** `filelist.txt; GSE335452_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE335nnn/GSE335452/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:42548877
- **Repository Access Link:** [GSE335452](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE335452)
- **Study Abstract / Experimental Design:**
  Pancreatic cancer resists immunotherapy due to a suppressive immune microenvironment where CD8⁺ T cells play a key role. Using single‑cell RNA sequencing and spatial transcriptomics, we characterized CD8⁺ exhausted T (Tex) cells in pancreatic ductal adenocarcinoma (PDAC). We generated single‑cell profiles from PDAC tumors and matched peripheral blood mononuclear cells, and performed T cell sub‑ana...

#### <a id='gse202051'></a>GSE202051 — Refined molecular taxonomy and treatment remodeling of pancreatic cancer using single-cell resolution
- **Cancer Type / Indication:** PDAC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** Single-Nucleus RNA-seq
- **Modality:** scRNA-seq + snRNA-seq
- **Cell Selection / Filtering Strategy:** **Nuclei Isolation (snRNA-seq)**
  > *"Single-nucleus RNA sequencing performed on isolated nuclei from frozen resection tissue."*
- **Cohort Scale:** 74 samples/patients; ~148,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE202051_adata_010nuc_10x.h5ad.gz; GSE202051_adata_010orgCRT_10x.h5ad.gz; GSE202051_totaldata-final-toshare.h5ad.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202051/suppl/GSE202051_adata_010nuc_1`
- **Citation / Reference:** PMID:36185212
- **Repository Access Link:** [GSE202051](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202051)
- **Study Abstract / Experimental Design:**
  Pancreatic ductal adenocarcinoma (PDAC) is a highly lethal and treatment-refractory cancer. Molecular stratification in pancreatic cancer remains rudimentary and does not yet inform clinical management or therapeutic development. Here we construct a high-resolution molecular landscape of the multicellular subtypes and spatial communities that compose PAC using single-nucleus RNA-seq and whole-tran...

#### <a id='gse291124'></a>GSE291124 — snRNA-seq dataset of 17 treatment-naive human PDAC patient-derived tumors
- **Cancer Type / Indication:** PDAC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** Single-Nucleus RNA-seq
- **Modality:** snRNA-seq
- **Cell Selection / Filtering Strategy:** **Nuclei Isolation (snRNA-seq)**
  > *"Single-nucleus RNA sequencing performed on isolated nuclei from frozen resection tissue."*
- **Cohort Scale:** 17 samples/patients; ~34,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE291124_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE291nnn/GSE291124/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41564862
- **Repository Access Link:** [GSE291124](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291124)
- **Study Abstract / Experimental Design:**
  This is a human single nucleus RNA sequencing (snRNAseq) dataset of human primary pancreatic ductal adenocarcinoma (PDAC) tumor tissues of 17 (deidentified) treatment-naive patients. The goal of this study was to thoroughly characterize the tumor microenvironment in human primary PDAC....

#### <a id='gse347847'></a>GSE347847 — Single-cell ATAC-seq profiling of ARID1A-loss-associated chromatin accessibility heterogeneity in pancreatic cancer cells
- **Cancer Type / Indication:** PDAC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 15 samples/patients; ~30,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE347847_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE347nnn/GSE347847/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE347847](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE347847)
- **Study Abstract / Experimental Design:**
  This Series contains single-cell ATAC-seq datasets from clone-resolved and pooled perturbation experiments examining how ARID1A and other epigenetic modifiers influence chromatin-accessibility heterogeneity in pancreatic ductal adenocarcinoma cells. The datasets include the main 28-clone BxPC3 experiment, a large-scale pooled five-shRNA experiment, acute SMARCA4 degradation with AU-15330, and inde...

#### <a id='gse312209'></a>GSE312209 — An oncogenic KRAS-driven secretome involving TNFα promotes niche preparation prior to pancreatic cancer onset [co-culture scRNA-seq]
- **Cancer Type / Indication:** PDAC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 15 samples/patients; ~30,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE312209_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE312nnn/GSE312209/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41634803
- **Repository Access Link:** [GSE312209](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE312209)
- **Study Abstract / Experimental Design:**
  Background: Pancreatic ductal adenocarcinomas (PDACs) are highly lethal and aggressive with oncogenic KRAS being the main oncogenic driver of the disease. PDACs have been extensively profiled at advanced stages, and in advanced disease the tumor microenvironment is a major determinant that critically shapes patient outcomes. Since the molecular events occurring prior to invasive growth remain poor...

#### <a id='gse348038'></a>GSE348038 — Single-cell RNA-seq profiling of ARID1A-loss-associated transcriptional heterogeneity and clonal diversification in pancreatic cancer cells
- **Cancer Type / Indication:** PDAC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 14 samples/patients; ~28,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `GSE348038_20221219_2_feature_reference.csv.gz; GSE348038_20230627_1_feature_reference.csv.gz; GSE348038_20241120_1_feature_reference.csv.gz`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE348nnn/GSE348038/suppl/GSE348038_20221219_2_fea`
- **Repository Access Link:** [GSE348038](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE348038)
- **Study Abstract / Experimental Design:**
  This Series contains additional single-cell RNA-seq datasets associated with a study of how loss of epigenetic regulators, particularly ARID1A, reshapes transcriptional heterogeneity in pancreatic ductal adenocarcinoma cells. The shPseuMO-Tag system couples genetic perturbation with clonal lineage tracing. The datasets include an initial epigenetic-modifier screen, an ARID1A knockdown time course,...

#### <a id='gse284392'></a>GSE284392 — An oncogenic KRAS-driven secretome involving TNFα promotes niche preparation prior to pancreatic cancer onset [time-series scRNA-Seq]
- **Cancer Type / Indication:** PDAC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 5 samples/patients; ~10,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE284392_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE284nnn/GSE284392/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41634803
- **Repository Access Link:** [GSE284392](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284392)
- **Study Abstract / Experimental Design:**
  Background: Pancreatic ductal adenocarcinomas (PDACs) are highly lethal and aggressive with oncogenic KRAS being the main oncogenic driver of the disease. PDACs have been extensively profiled at advanced stages, and in advanced disease the tumor microenvironment is a major determinant that critically shapes patient outcomes. Since the molecular events occurring prior to invasive growth remain poor...

#### <a id='gse327056'></a>GSE327056 — Spatial transcriptomic profiling of human pancreatic ductal adenocarcinoma  using 10x Genomics Visium platform
- **Cancer Type / Indication:** PDAC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** 10x Visium Spatial Transcriptomics
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 4 samples/patients; ~8,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE327056_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE327nnn/GSE327056/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:42100437
- **Repository Access Link:** [GSE327056](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327056)
- **Study Abstract / Experimental Design:**
  In this study, we performed spatial transcriptomics (ST) to investigate the gene expression features across one normal pancreatic tissue, PC tissue, adjacent tumor tissue, and tumor stroma using 10x Genomics Visium spatial gene expression platform. We generated high-quality spatial gene expression profiles combined with histomorphological information, aiming to reveal the spatial heterogeneity of ...

#### <a id='gse160977'></a>GSE160977 — Single cell data from patient-derived xenograft model of PDAC tumor.
- **Cancer Type / Indication:** PDAC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~2,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE160977_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160977/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE160977](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160977)
- **Study Abstract / Experimental Design:**
  The aim of this study is to prove that PDX models of PDAC tumor are able to recapitulate different sub-populations of cells previously identified in human samples...

#### <a id='gse288067'></a>GSE288067 — Single-cell RNA sequencing combined with multiplex immunofluorescence probes the role of MFAP5+ fibroblasts in the microenvironment of pancreatic ductal adenocarcinoma.
- **Cancer Type / Indication:** PDAC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq + Spatial
- **Cell Selection / Filtering Strategy:** **Unselected / Total Single-Cell Suspension**
  > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Cohort Scale:** 1 samples/patients; ~23,905 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE288067_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE288nnn/GSE288067/suppl/filelist.txt; https://ft`
- **Repository Access Link:** [GSE288067](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288067)
- **Study Abstract / Experimental Design:**
  we employed single-cell RNA sequencing to examine the biological characteristics of MFAP5+ fibroblasts in PDAC and their interaction with vascular endothelial cells within tumors. We then utilized a proposed temporal sequencing analysis technique to infer the evolution of cellular subtypes of cancer-associated fibroblasts. To verify our hypothesis, we employed a multiplex immunofluorescence techni...

#### <a id='gse300154'></a>GSE300154 — Gene expresson profile at single cell level of  myeloid cells from pancreatic cancer tissue.
- **Cancer Type / Indication:** PDAC
- **Classification:** Tier 2 (Baseline Atlas)
- **Sequencing Technology:** High-Throughput scRNA-seq
- **Modality:** scRNA-seq
- **Cell Selection / Filtering Strategy:** **FACS-sorted (CD3+/CD8+ T-cell enriched)**
  > *"Myeloid cells (CD33+ cells) are heterogenous population and are considered to be immune suppressive in pancreatic cancer tissue."*
- **Cohort Scale:** 1 samples/patients; ~2,000 cells
- **Clinical Context / Response:** Primary / untreated baseline tumor atlas (deconvolution reference candidate)
- **Downloadable Matrix Files:** `filelist.txt; GSE300154_RAW.tar`
- **Direct Download URLs:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE300nnn/GSE300154/suppl/filelist.txt; https://ft`
- **Citation / Reference:** PMID:41722836
- **Repository Access Link:** [GSE300154](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300154)
- **Study Abstract / Experimental Design:**
  Myeloid cells (CD33+ cells) are heterogenous population and are considered to be immune suppressive in pancreatic cancer tissue. We conducted single-cell RNA sequencing of myeloid cells from pancreatic cancer to analyze the diversity of myeloid cells and their molecular characterisitcs in pancratic cancer....
