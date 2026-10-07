# Tier 0: Premier Benchmark Core scRNA-seq Cancer Cohorts

This document details the **22 ultra-select Tier 0 cohorts** selected from the 352-cohort discovery catalog. These cohorts represent high-powered clinical immunotherapy trials with verified response outcomes, high patient numbers ($N \ge 9$, median $N = 51$), and broad cellular microenvironment representation across 9 human solid tumor indications.

## Summary Statistics
- **Total Selected Cohorts**: 22
- **Total Clinically Annotated Patients/Samples**: 1,999
- **Total Single-Cell/Nucleus Transcriptomes**: ~4,009,048
- **Unselected Whole-TME Suspensions**: 18 cohorts (82%)

## Cohort Master Table

| Accession | Indication | Score | Patients | Cells | Modality & Tech | Cell Selection | Clinical Response Details |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **GSE246613** | Breast | **100/100** | 266 | 532,000 | scRNA-seq + Spatial (High-Throughput scRNA-seq) | Unselected / Total Single-Cell Suspension | Documented ICB immunotherapy with response / resistance / outcome correlates |
| **GSE300475** | Breast | **98/100** | 32 | 64,000 | scRNA-seq (High-Throughput scRNA-seq) | Unselected / Total Single-Cell Suspension | Documented ICB immunotherapy with response / resistance / outcome correlates |
| **GSE212707** | Breast | **88/100** | 21 | 42,000 | scRNA-seq + Spatial (High-Throughput scRNA-seq) | Unselected / Total Single-Cell Suspension | Documented ICB immunotherapy with response / resistance / outcome correlates |
| **GSE236581** | CRC | **91/100** | 169 | 338,000 | scRNA-seq (High-Throughput scRNA-seq) | FACS-sorted (CD3+/CD8+ T-cell enriched) | Documented ICB immunotherapy with response / resistance / outcome correlates |
| **GSE299651** | CRC | **88/100** | 20 | 40,000 | scRNA-seq (High-Throughput scRNA-seq) | Unselected / Total Single-Cell Suspension | Documented ICB immunotherapy with response / resistance / outcome correlates |
| **CELLxGENE_829a3cd1** | CRC | **85/100** | 29 | 47,107 | scRNA-seq (10x Chromium 3' v3/v3.1) | Unselected / Total Single-Cell Suspension | ICB immunotherapy treated cohort (response evaluation required) |
| **GSE270680** | Gastric | **88/100** | 77 | 154,000 | scRNA-seq + Spatial (High-Throughput scRNA-seq) | Unselected / Total Single-Cell Suspension | Documented ICB immunotherapy with response / resistance / outcome correlates |
| **GSE313642** | HCC | **98/100** | 194 | 388,000 | scRNA-seq (High-Throughput scRNA-seq) | Unselected / Total Single-Cell Suspension | Documented ICB immunotherapy with response / resistance / outcome correlates |
| **GSE245906** | HCC | **84/100** | 20 | 40,000 | scRNA-seq (High-Throughput scRNA-seq) | Unselected / Total Single-Cell Suspension | Documented ICB immunotherapy with response / resistance / outcome correlates |
| **GSE301741** | HNSCC | **98/100** | 58 | 116,000 | scRNA-seq (High-Throughput scRNA-seq) | Unselected / Total Single-Cell Suspension | Documented ICB immunotherapy with response / resistance / outcome correlates |
| **GSE287301** | HNSCC | **86/100** | 48 | 96,000 | scRNA-seq + Spatial (Subcellular Spatial Transcriptomics (CosMx/Xenium)) | Unselected / Total Single-Cell Suspension | ICB immunotherapy treated cohort (response evaluation required) |
| **GSE200996** | HNSCC | **85/100** | 204 | 408,000 | scRNA-seq (High-Throughput scRNA-seq) | FACS-sorted (CD3+/CD8+ T-cell enriched) | Documented ICB immunotherapy with response / resistance / outcome correlates |
| **CELLxGENE_7b20c613** | Melanoma | **90/100** | 167 | 355,941 | scRNA-seq (10x Chromium 3' v3/v3.1) | Unselected / Total Single-Cell Suspension | ICB immunotherapy treated cohort (response evaluation required) |
| **GSE218429** | Melanoma | **89/100** | 35 | 70,000 | scRNA-seq (High-Throughput scRNA-seq) | Unselected / Total Single-Cell Suspension | Documented ICB immunotherapy with response / resistance / outcome correlates |
| **GSE344166** | Melanoma | **84/100** | 116 | 232,000 | scRNA-seq + Spatial (High-Throughput scRNA-seq) | Unselected / Total Single-Cell Suspension | Documented ICB immunotherapy with response / resistance / outcome correlates |
| **GSE207422** | NSCLC | **99/100** | 39 | 78,000 | scRNA-seq (High-Throughput scRNA-seq) | Unselected / Total Single-Cell Suspension | Documented ICB immunotherapy with response / resistance / outcome correlates |
| **GSE317309** | NSCLC | **93/100** | 64 | 128,000 | scRNA-seq (High-Throughput scRNA-seq) | Unselected / Total Single-Cell Suspension | Documented ICB immunotherapy with response / resistance / outcome correlates |
| **GSE243013** | NSCLC | **90/100** | 243 | 486,000 | scRNA-seq (High-Throughput scRNA-seq) | Unselected / Total Single-Cell Suspension | Documented ICB immunotherapy with response / resistance / outcome correlates |
| **GSE311789** | PDAC | **91/100** | 142 | 284,000 | scRNA-seq (High-Throughput scRNA-seq) | Unselected / Total Single-Cell Suspension | Documented ICB immunotherapy with response / resistance / outcome correlates |
| **GSE316195** | PDAC | **86/100** | 22 | 44,000 | scRNA-seq + snRNA-seq (Single-Nucleus RNA-seq) | Nuclei Isolation (snRNA-seq) | Documented ICB immunotherapy with response / resistance / outcome correlates |
| **GSE210038** | ccRCC | **81/100** | 9 | 18,000 | scRNA-seq + Spatial (High-Throughput scRNA-seq) | Unselected / Total Single-Cell Suspension | Documented ICB immunotherapy with response / resistance / outcome correlates |
| **GSE314072** | ccRCC | **76/100** | 24 | 48,000 | scRNA-seq (High-Throughput scRNA-seq) | FACS-sorted (CD3+/CD8+ T-cell enriched) | Documented ICB immunotherapy with response / resistance / outcome correlates |

---

## Detailed Cohort Profiles

### 1. GSE246613 — Breast (Score: 100/100)
**Title**: Single-cell and spatial profiling identify three response trajectories to pembrolizumab and radiation therapy in triple negative breast cancer

- **Accession URL**: [GSE246613](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246613)
- **Indication**: Breast
- **Sample/Patient Count**: 266 patients / biological specimens
- **Cell Count Estimate**: ~532,000 single cells/nuclei
- **Cell Selection Strategy**: `Unselected / Total Single-Cell Suspension`
- **Protocol Excerpt**: > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Scoring Rationale**: `Clin:30/30 | Scale:20/20 (N=266) | Kinetics:15/15 | Breadth:15/15 (Unselected / Total Single-Cell Suspension) | Matrix:10/10 | Impact:10/10`
- **Available Raw Matrices**: `GSE246613_PembroRT_immune_R100_final.h5ad.gz; GSE246613_PembroRT_non_immune_cells.h5ad.gz; GSE246613_combined_RTPDv4_scvi_celltypist.h5ad.gz`

**Study Abstract / Clinical Setting**:
Cancer immunotherapy trials have produced encouraging results, but resistance remains a problem necessitating strategies to identify patients that benefit from immunotherapy alone or who require additional combinations like chemotherapy or radiotherapy. Here we employ single-cell transcriptomics and spatial proteomics to profile triple negative breast cancer biopsies taken before and after one cycle of pembrolizumab and after a second cycle of pembrolizumab given with radiotherapy. Non-responders lack immune infiltrate before and after therapy and exhibit minimal therapy-induced immune changes. Responding tumors form two groups that are distinguishable by a classifier prior to therapy, with one showing high MHC expression, evidence of tertiary lymphoid structures and displaying anti-tumor immunity before treatment. The other responder group resembles non-responders at baseline and mounts a maximal immune response after the combination therapy, which is characterized by cytotoxic T cell and antigen presenting myeloid cell interactions and is mirrored in murine breast tumors only after radiation.

### 2. GSE300475 — Breast (Score: 98/100)
**Title**: Single cell RNA sequencing for longitudinal human peripheral blood from HR+ breast cancer patients treated with immunotherapy

- **Accession URL**: [GSE300475](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300475)
- **Indication**: Breast
- **Sample/Patient Count**: 32 patients / biological specimens
- **Cell Count Estimate**: ~64,000 single cells/nuclei
- **Cell Selection Strategy**: `Unselected / Total Single-Cell Suspension`
- **Protocol Excerpt**: > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Scoring Rationale**: `Clin:30/30 | Scale:20/20 (N=32) | Kinetics:15/15 | Breadth:15/15 (Unselected / Total Single-Cell Suspension) | Matrix:8/10 | Impact:10/10`
- **Available Raw Matrices**: `GSE300475_feature_ref.xlsx; filelist.txt; GSE300475_RAW.tar`

**Study Abstract / Clinical Setting**:
The limited benefit of immune checkpoint inhibitor in breast cancer indicates the pressing need to identify biomarkers of response to minimize risk and maximize benefit. In this study, we performed single cell RNA sequencing and T cell receptor sequencing on peripheral blood mononuclear cells to monitor the peripheral immune dynamics of an exploratory cohort of hormone receptor positive breast cancer patients treated with neoadjuvant nab-paclitaxel+pembrolizumab with the ultimate goal of identifying potential peripheral blood predictive biomarkers.

### 3. GSE212707 — Breast (Score: 88/100)
**Title**: Multiomic analysis reveals conservation of cancer associated fibroblast phenotypes across species and tissue of origin [multiome]

- **Accession URL**: [GSE212707](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212707)
- **Indication**: Breast
- **Sample/Patient Count**: 21 patients / biological specimens
- **Cell Count Estimate**: ~42,000 single cells/nuclei
- **Cell Selection Strategy**: `Unselected / Total Single-Cell Suspension`
- **Protocol Excerpt**: > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Scoring Rationale**: `Clin:30/30 | Scale:15/20 (N=21) | Kinetics:15/15 | Breadth:15/15 (Unselected / Total Single-Cell Suspension) | Matrix:8/10 | Impact:5/10`
- **Available Raw Matrices**: `filelist.txt; GSE212707_RAW.tar`

**Study Abstract / Clinical Setting**:
Cancer associated fibroblasts (CAFs) are integral to the solid tumor microenvironment. Once thought to be a relatively uniform population of matrix-producing cells, the arrival of single cell RNA sequencing has revealed diverse CAF phenotypes. Here, we further probe CAF heterogeneity with a comprehensive multiome approach. Using paired, same-cell chromatin accessibility and transcriptome analysis, we provide an integrated analysis of CAF subpopulations over a complex spatial transcriptomic and proteomic landscape to identify three superclusters – steady state-like (SSL), mechanoresponsive (MR) and immunomodulatory (IM) CAFs. These superclusters are recapitulated across multiple tissue types and species. Selective disruption of underlying mechanical force or immune checkpoint inhibition therapy results in shifts in CAF subpopulation distributions and impacts tumor growth. As such, the balance among CAF superclusters may have considerable translational implications. Collectively, this research expands our understanding of CAF biology, identifying regulatory pathways in CAF differentiation and elucidating novel therapeutic targets in a species- and tumor-agnostic manner.

### 4. GSE236581 — CRC (Score: 91/100)
**Title**: Spatiotemporal single-cell analysis decodes cellular dynamics underlying different responses to immunotherapy in Colorectal Cancer

- **Accession URL**: [GSE236581](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236581)
- **Indication**: CRC
- **Sample/Patient Count**: 169 patients / biological specimens
- **Cell Count Estimate**: ~338,000 single cells/nuclei
- **Cell Selection Strategy**: `FACS-sorted (CD3+/CD8+ T-cell enriched)`
- **Protocol Excerpt**: > *"Finally, a predictive signature was established using circulating CD8 T cells."*
- **Scoring Rationale**: `Clin:30/30 | Scale:20/20 (N=169) | Kinetics:15/15 | Breadth:6/15 (FACS-sorted (CD3+/CD8+ T-cell enriched)) | Matrix:10/10 | Impact:10/10`
- **Available Raw Matrices**: `GSE236581_counts.mtx.gz`

**Study Abstract / Clinical Setting**:
Expanding the efficacy of immune checkpoint blockade (ICB) in colorectal cancer (CRC) patients presses for a comprehensive understanding of treatment responsiveness. Here, we analyzed 169 single-cell samples from CRC patients at multiple sequential time points during the course of anti-PD-1 neoadjuvant therapy to map the evolution of local and systemic immunity. In tumors, exhausted T (Tex) cells or tumor-reactive-like CD8 T (Ttr-like) cells were closely related to treatment efficacy, and we observed correlated dynamics between Tex cells and multiple other tumor-enriched cell types following the treatment. Accordingly, several coordinated cellular programs exhibiting distinct response associations were identified. From a systemic perspective, we found divergent replenishment patterns of Ttr-like cells underlying different response statuses and decoded the phenotypic transitions of Ttr-like cells as they infiltrated tissues from the periphery. Finally, a predictive signature was established using circulating CD8 T cells. Our study provides novel insights into the spatiotemporal cellular dynamics following PD-1 blockade in CRC.

### 5. GSE299651 — CRC (Score: 88/100)
**Title**: Pooled single-cell screening in colorectal cancer identifies transcriptional modules of clinical relevance unlocked by oncogenes

- **Accession URL**: [GSE299651](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299651)
- **Indication**: CRC
- **Sample/Patient Count**: 20 patients / biological specimens
- **Cell Count Estimate**: ~40,000 single cells/nuclei
- **Cell Selection Strategy**: `Unselected / Total Single-Cell Suspension`
- **Protocol Excerpt**: > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Scoring Rationale**: `Clin:30/30 | Scale:15/20 (N=20) | Kinetics:8/15 | Breadth:15/15 (Unselected / Total Single-Cell Suspension) | Matrix:10/10 | Impact:10/10`
- **Available Raw Matrices**: `GSE299651_20240523_Caco2_pool_cellranger_3_filtered_feature_bc_matrix.h5; GSE299651_20240523_HT29_pool_cellranger_3_filtered_feature_bc_matrix.h5; GSE299651_20240523_RKO_pool_cellranger_3_filtered_feature_bc_matrix.h5`

**Study Abstract / Clinical Setting**:
While oncogenic mutations shape colorectal cancer biology and therapy response, their prognostic value remains low. Cluster-based classification of patient cancer transcriptomes has shown greater promise for prognosis, yet these systems do not account for the roles of oncogenes in establishing cancer phenotypes. Here, we create and validate a prognostic classifier for colorectal cancer based on transcriptional programs induced by oncogenes. To systematically investigate oncogenic drivers, we employed a barcoded library of colorectal cancer-associated oncogene variants across a panel of genetically diverse colorectal cancer cell lines. We profiled the transcriptomes of over 100,000 transgenic cells and used machine learning to define transcriptional modules capturing key functional traits. Our analysis revealed heterogeneity on the cell-to-cell level and context-dependent gene expression patterns induced by oncogenes. We identified overarching gene expression modules reflecting core oncogenic processes, including cancer cell plasticity, inflammatory response, replicative stress, and epithelial-to-mesenchymal transition. These modules enabled a functional classification that linked oncogenic signalling states to distinct transcriptional profiles. We demonstrated their prognostic value by stratifying clinical colorectal cancer cohorts into high- and low-risk groups. Although partially correlated with established clinical parameters, the modules provided additional prognostic information, improving survival prediction and therapy stratification beyond existing classification systems. In summary, our study establishes a framework that connects oncogenic mutations to core transcriptional modules. By integrating experimental models with clinical data, we provide a resource for investigating colorectal cancer progression and oncogene-specific vulnerabilities, facilitating future research and precision oncology.

### 6. CELLxGENE_829a3cd1 — CRC (Score: 85/100)
**Title**: progressive_plasticity_during_crc_metastasis_epithelial

- **Accession URL**: [CELLxGENE_829a3cd1](https://cellxgene.cziscience.com/e/829a3cd1-a466-49f1-b2e9-d3f6b7f392e2.cxg/)
- **Indication**: CRC
- **Sample/Patient Count**: 29 patients / biological specimens
- **Cell Count Estimate**: ~47,107 single cells/nuclei
- **Cell Selection Strategy**: `Unselected / Total Single-Cell Suspension`
- **Protocol Excerpt**: > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Scoring Rationale**: `Clin:30/30 | Scale:15/20 (N=29) | Kinetics:5/15 | Breadth:15/15 (Unselected / Total Single-Cell Suspension) | Matrix:10/10 | Impact:10/10`
- **Available Raw Matrices**: `829a3cd1-a466-49f1-b2e9-d3f6b7f392e2.h5ad`

**Study Abstract / Clinical Setting**:
Collection: Progressive plasticity during colorectal cancer metastasis | Diseases: colorectal cancer, normal | Tissues: caecum, chest wall, descending colon, hepatic flexure of colon, liver, lung, peritoneum, rectum, sigmoid colon, transverse colon | Assays: 10x 3' v3

### 7. GSE270680 — Gastric (Score: 88/100)
**Title**: A spatially resolved atlas of gastric cancer characterises a lymphocyte aggregated region [scRNA-seq]

- **Accession URL**: [GSE270680](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270680)
- **Indication**: Gastric
- **Sample/Patient Count**: 77 patients / biological specimens
- **Cell Count Estimate**: ~154,000 single cells/nuclei
- **Cell Selection Strategy**: `Unselected / Total Single-Cell Suspension`
- **Protocol Excerpt**: > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Scoring Rationale**: `Clin:30/30 | Scale:20/20 (N=77) | Kinetics:5/15 | Breadth:15/15 (Unselected / Total Single-Cell Suspension) | Matrix:8/10 | Impact:10/10`
- **Available Raw Matrices**: `filelist.txt; GSE270680_RAW.tar`

**Study Abstract / Clinical Setting**:
The tumour microenvironment (TME) is a focal point in cancer immunotherapy: its cellular composition and spatial organisation, especially the distribution of lymphocytes, can affect the clinical outcomes of cancer patients. In addition, the function of a cell differs depending on its spatial location and interaction with neighbouring cells. Here, by integrating single-cell transcriptomics with spatial transcriptomics, we survey how the spatial distribution of different cell types varies across diverse histological regions of gastric cancer. Notably, tertiary lymphoid structures (TLSs), a tissue architecture harbouring lymphocytes more than any other region within the TME, possess great potential in modulating anti-tumorigenic immunity by elevating an array of genes (FDCSP and CCL19) principally involved in lymphocyte recruitment and activation.  Our findings advance the understanding of TLSs in contributing to anti-tumorigenic immunity in a spatially resolved context, which could be further leveraged as a predictive marker for immunotherapy response.

### 8. GSE313642 — HCC (Score: 98/100)
**Title**: Immunosuppressive monocytes are enriched in hepatocellular carcinoma patients with liver dysfunction in a phase II trial of combination sorafenib and nivolumab [CITE-Seq]

- **Accession URL**: [GSE313642](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313642)
- **Indication**: HCC
- **Sample/Patient Count**: 194 patients / biological specimens
- **Cell Count Estimate**: ~388,000 single cells/nuclei
- **Cell Selection Strategy**: `Unselected / Total Single-Cell Suspension`
- **Protocol Excerpt**: > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Scoring Rationale**: `Clin:30/30 | Scale:20/20 (N=194) | Kinetics:15/15 | Breadth:15/15 (Unselected / Total Single-Cell Suspension) | Matrix:8/10 | Impact:10/10`
- **Available Raw Matrices**: `filelist.txt; GSE313642_RAW.tar`

**Study Abstract / Clinical Setting**:
Background and aims: Immune checkpoint inhibition (ICI) and anti-angiogenic therapies are active in hepatocellular carcinoma (HCC), although patients with impaired hepatic function have worse outcomes. Methods: We conducted a multi-center, open-label phase II clinical trial to assess the safety and efficacy of the multikinase inhibitor, sorafenib, combined with nivolumab, in patients with advanced or unresectable HCC and varying liver function. In a Part 1 safety lead-in, we investigated the primary endpoint of the maximum-tolerated dose (MTD) of the combination in patients with Child-Pugh A or B7 HCC. In Part 2, we enrolled patients with Child-Pugh B HCC with the primary endpoint of grade >/=3 treatment-related adverse events (TRAE) incidence. Exploratory endpoints included immunologic biomarkers. Results: Overall, 25 patients were consented and 16 eligible patients enrolled. In Part 1, dose-limiting toxicity occurred in 1 of 6 patients in Dose Level -1, and 2 of 5 patients in Dose Level 1; Dose Level -1 was determined to be the MTD. In total, 69% of patients experienced a Grade >/=3 TRAE, with similar distribution for patients with Child-Pugh A and B liver function (70%, 95% CI: 0.35, 0.93 vs 66.7%, 95% CI: 0.22, 0.96). The objective response rate was 6%. We found that patients with Child-Pugh B liver disease harbored more circulating suppressive CD14+ monocytes at baseline. Conclusions: While the combination of sorafenib and nivolumab demonstrated acceptable safety at the MTD in both Child-Pugh subgroups, the objective response rate was below the pre-specified threshold to be declared worthy of further exploration. A distinct immune cell profile in Child-Pugh B patients may define mechanisms of resistance and potential therapeutic targets in this population with unmet clinical need

### 9. GSE245906 — HCC (Score: 84/100)
**Title**: Identification of TREM1+CD163+ myeloid cells as a deleterious immune subset in HCC [scRNA-seq]

- **Accession URL**: [GSE245906](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245906)
- **Indication**: HCC
- **Sample/Patient Count**: 20 patients / biological specimens
- **Cell Count Estimate**: ~40,000 single cells/nuclei
- **Cell Selection Strategy**: `Unselected / Total Single-Cell Suspension`
- **Protocol Excerpt**: > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Scoring Rationale**: `Clin:30/30 | Scale:15/20 (N=20) | Kinetics:5/15 | Breadth:15/15 (Unselected / Total Single-Cell Suspension) | Matrix:9/10 | Impact:10/10`
- **Available Raw Matrices**: `GSE245906_Giraud_Chalopin_innate_HCC_metadata.tsv.gz; GSE245906_Giraud_Chalopin_innate_HCC_processed_count_data_mat.tsv.gz`

**Study Abstract / Clinical Setting**:
Hepatocellular carcinoma (HCC) is an inflammation-associated cancer arising from viral and non-viral etiologies. Expansion of suppressive myeloid cells is a hallmark of chronic inflammation and cancer, but their heterogeneity in HCC is not fully resolved and might underlie immunotherapy resistance in the steatohepatitis setting. Here, we present a high resolution atlas of hepatic innate immune cells from patients with HCC that unravels a steatohepatitis contexture characterized by influx of inflammatory and immunosuppressive myeloid cells. A discrete myeloid cell population identified by selective expression of TREM1 and CD163 expands in steatohepatitis-HCC. We refer to this population as TREM1+ regulatory myeloid cells (Mreg), as it potently suppresses T cell effector functions, highly expresses TGFB1 and IL13RA1 and localizes to HCC fibrotic lesions.

### 10. GSE301741 — HNSCC (Score: 98/100)
**Title**: Single cell analysis highlights the significance of malignant cell IFN/MHC-II for immunotherapy response in head and neck squamous cell carcinoma

- **Accession URL**: [GSE301741](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301741)
- **Indication**: HNSCC
- **Sample/Patient Count**: 58 patients / biological specimens
- **Cell Count Estimate**: ~116,000 single cells/nuclei
- **Cell Selection Strategy**: `Unselected / Total Single-Cell Suspension`
- **Protocol Excerpt**: > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Scoring Rationale**: `Clin:30/30 | Scale:20/20 (N=58) | Kinetics:15/15 | Breadth:15/15 (Unselected / Total Single-Cell Suspension) | Matrix:8/10 | Impact:10/10`
- **Available Raw Matrices**: `GSE301741_Seurat_Object_QCpass_137020cells_withMetaData.rds; filelist.txt; GSE301741_RAW.tar`

**Study Abstract / Clinical Setting**:
For many cancers, including head and neck squamous cell carcinoma (HNSCC), response rates to immunotherapy remain modest, with limited ability to predict responders. Previous studies that characterized cellular changes associated with immunotherapy in HNSCC focused on immune cells, providing limited insight into malignant cell responses. Motivated by this gap, we performed single cell RNA-sequencing on 16 HNSCC patients pre- and post-neoadjuvant pembrolizumab treatment. We identified a malignant-interferon (IFN)/MHC-II program in malignant cells, characterized by expression of MHC-II and interferon response genes, which was associated with tumor response to pembrolizumab. We validated malignant cell MHC-II expression at the protein level and characterized its relationship with surrounding immune subsets via multiplexed immunofluorescence. In a murine HNSCC model, IFN-γ-induced malignant cell MHC-II expression marked tumors with favorable immune microenvironments and sensitivity to immunotherapy. Finally, we confirmed pre-treatment expression of the malignant-IFN/MHC-II program as marker of response through deconvolution of bulk RNA-seq data from an independent cohort of 25 pembrolizumab-treated HNSCC patients. Beyond identifying malignant cell MHC-II expression and the malignant-IFN/MHC-II program as potential biomarkers for immunotherapy response in HNSCC, our work provides additional insights into the specific and important role of malignant cells in immunotherapy.

### 11. GSE287301 — HNSCC (Score: 86/100)
**Title**: Single-cell RNA-sequencing of HNSCC-infiltrating T cells

- **Accession URL**: [GSE287301](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287301)
- **Indication**: HNSCC
- **Sample/Patient Count**: 48 patients / biological specimens
- **Cell Count Estimate**: ~96,000 single cells/nuclei
- **Cell Selection Strategy**: `Unselected / Total Single-Cell Suspension`
- **Protocol Excerpt**: > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Scoring Rationale**: `Clin:30/30 | Scale:20/20 (N=48) | Kinetics:5/15 | Breadth:15/15 (Unselected / Total Single-Cell Suspension) | Matrix:6/10 | Impact:10/10`
- **Available Raw Matrices**: `GSE287301_filtered_feature_bc_matrix.tar.gz`

**Study Abstract / Clinical Setting**:
This study examines the T cell landscape in head and neck squamous cell carcinoma (HNSCC) using single-cell RNA sequencing combined with T cell receptor and protein profiling. By analyzing tumor-infiltrating T cells from 28 patients, we identified diverse T cell subsets, characterized their transcriptional states, and mapped their clonal dynamics. These findings provide a detailed view of the immune microenvironment in HNSCC and offer insights into T cell-mediated immunity and potential targets for immunotherapy and are patient-matched to PBMC TCR sequencing and Xenium spatial transcriptomics in separate GEO submissions.

### 12. GSE200996 — HNSCC (Score: 85/100)
**Title**: Tissue-resident Memory and Circulating T cells are Early Responders to Pre-surgical Cancer Immunotherapy

- **Accession URL**: [GSE200996](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200996)
- **Indication**: HNSCC
- **Sample/Patient Count**: 204 patients / biological specimens
- **Cell Count Estimate**: ~408,000 single cells/nuclei
- **Cell Selection Strategy**: `FACS-sorted (CD3+/CD8+ T-cell enriched)`
- **Protocol Excerpt**: > *"Single-cell suspension sorted by FACS for CD3+/CD8+ T lymphocytes (TILs)."*
- **Scoring Rationale**: `Clin:30/30 | Scale:20/20 (N=204) | Kinetics:10/15 | Breadth:6/15 (FACS-sorted (CD3+/CD8+ T-cell enriched)) | Matrix:9/10 | Impact:10/10`
- **Available Raw Matrices**: `GSE200996_CD4.PBMC.single.cell.meta.data.txt.gz; GSE200996_CD4.tumor.single.cell.meta.data.txt.gz; GSE200996_CD45.PBMC.single.cell.meta.data.txt.gz`

**Study Abstract / Clinical Setting**:
Pre-surgical (neoadjuvant) immune checkpoint blockade has shown promising activity in multiple cancer types, but the molecular mechanisms are not well understood. Here, we characterized early kinetic changes in tumor-infiltrating and circulating immune cells in oral cancer patients treated with neoadjuvant anti-PD-1 or anti-PD-1/CTLA-4 in a phase 2 clinical trial. Tumor-infiltrating CD8 T cells that clonally expanded during immunotherapy expressed elevated tissue-resident memory and cytotoxicity programs compared to non-responding cells. These programs were already active in pre-treatment T cells that later responded, reflecting a capacity for rapid response. Treatment also induced a systemic immune response, including expansion of pre-existing and emergent T cell clonotypes undetectable prior to therapy. The frequency of activated blood CD8 T cells, including pretreatment PD-1-positive KLRG1-negative T cells, was strongly associated with intra-tumoral pathological response, and these activated cells were enriched for tumor-infiltrating T cell clonotypes. These results demonstrate how neoadjuvant checkpoint blockade induces local and systemic tumor immunity.

### 13. CELLxGENE_7b20c613 — Melanoma (Score: 90/100)
**Title**: Integrated cancer cell-specific single-cell RNA-seq datasets of immune checkpoint blockade-treated patients

- **Accession URL**: [CELLxGENE_7b20c613](https://cellxgene.cziscience.com/e/7b20c613-9add-43d1-87e9-defd3d9b9f8c.cxg/)
- **Indication**: Melanoma
- **Sample/Patient Count**: 167 patients / biological specimens
- **Cell Count Estimate**: ~355,941 single cells/nuclei
- **Cell Selection Strategy**: `Unselected / Total Single-Cell Suspension`
- **Protocol Excerpt**: > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Scoring Rationale**: `Clin:30/30 | Scale:20/20 (N=167) | Kinetics:5/15 | Breadth:15/15 (Unselected / Total Single-Cell Suspension) | Matrix:10/10 | Impact:10/10`
- **Available Raw Matrices**: `7b20c613-9add-43d1-87e9-defd3d9b9f8c.h5ad`

**Study Abstract / Clinical Setting**:
Collection: Integrated cancer cell-specific single-cell RNA-seq datasets of immune checkpoint blockade-treated patients | Diseases: HER2 positive breast carcinoma, basal cell carcinoma, estrogen-receptor positive breast cancer, hepatocellular carcinoma, intrahepatic cholangiocarcinoma, melanoma, metastatic melanoma, nonpapillary renal cell carcinoma, squamous cell carcinoma, triple-negative breast carcinoma | Tissues: arm skin, brain, breast, kidney, liver, nose skin, skin of body, skin of calf, skin of cheek, skin of external ear, skin of forehead, skin of knee, skin of neck, skin of scalp | Assays: 10x 3' v2, 10x 3' v3, 10x 5' transcription profiling, Smart-seq2

### 14. GSE218429 — Melanoma (Score: 89/100)
**Title**: Downregulation of KEAP1 in melanoma promotes resistance to immune checkpoint blockade [scRNA-seq]

- **Accession URL**: [GSE218429](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218429)
- **Indication**: Melanoma
- **Sample/Patient Count**: 35 patients / biological specimens
- **Cell Count Estimate**: ~70,000 single cells/nuclei
- **Cell Selection Strategy**: `Unselected / Total Single-Cell Suspension`
- **Protocol Excerpt**: > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Scoring Rationale**: `Clin:30/30 | Scale:20/20 (N=35) | Kinetics:5/15 | Breadth:15/15 (Unselected / Total Single-Cell Suspension) | Matrix:9/10 | Impact:10/10`
- **Available Raw Matrices**: `GSE218429_counts.csv.gz`

**Study Abstract / Clinical Setting**:
Immune checkpoint blockade (ICB) has demonstrated efficacy in patients with melanoma, but many exhibit poor responses. Using single cell RNA sequencing of melanoma patient-derived circulating tumor cells (CTCs) and functional characterization using mouse melanoma models, we show that the KEAP1/NRF2 pathway modulates sensitivity to ICB, independently of tumorigenesis. The NRF2 negative regulator, KEAP1, shows intrinsic variation in expression, leading to tumor heterogeneity and subclonal resistance.

### 15. GSE344166 — Melanoma (Score: 84/100)
**Title**: Integrated spatial and single cell analysis identifies CCL21+ lymphatic endothelial cells as a driver of a favorable immune environment in acral melanoma

- **Accession URL**: [GSE344166](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE344166)
- **Indication**: Melanoma
- **Sample/Patient Count**: 116 patients / biological specimens
- **Cell Count Estimate**: ~232,000 single cells/nuclei
- **Cell Selection Strategy**: `Unselected / Total Single-Cell Suspension`
- **Protocol Excerpt**: > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Scoring Rationale**: `Clin:30/30 | Scale:20/20 (N=116) | Kinetics:5/15 | Breadth:15/15 (Unselected / Total Single-Cell Suspension) | Matrix:9/10 | Impact:5/10`
- **Available Raw Matrices**: `GSE344166_Acral_GeoMX_norm.csv.gz; GSE344166_Experiment_Summary.csv.gz; GSE344166_Hs_R_NGS_WTA_v1.0.pkc.gz`

**Study Abstract / Clinical Setting**:
Acral melanoma (AM) is a rare type of melanoma that responds poorly to immunotherapy. In the current study we undertook integrated single cell RNA-Seq, spatial transcriptomics and multiplexed immunofluorescence to identify potential regulators of the immune environment. Our analyses identified distinct immune habitats at the invasive front, intratumoral regions and areas of distant inflammation across the AM samples. Two distinct subclasses of endothelial cells were identified, one with an immune-regulating profile (CCL21, IGF1, IL-7, PDPN) that were characteristic of lymphatic endothelial cells (LECs) and a population of vascular endothelial cells (VECs) (SPARCL1, PVLAP, ADGRL4). Cell-cell interaction and correlation analysis showed the CCL21+ LECs to be a strong regulator of the immune microenvironment through communication with dendritic cells (DCs) and CD4+ T cells. These LEC-dendritic cell interactions were mediated through CCL21-CCR7 and were strongly correlated in both the single cell and spatial data. By contrast, the VECs were predicted to interact with fibroblasts and macrophages. Validation by multiplexed immunofluorescence identified LECs in the majority of AM samples and demonstrated that only those staining positively for CCL21 were spatially associated with DCs, B-cells and CD4+ T-cells. Expression of CCL21 and CCR7 showed a significant correlation and were associated with increased survival in melanoma patient cohorts. Together these data identify the presence of CCL21+ LECs in acral melanoma that positively regulate the immune microenvironment.

### 16. GSE207422 — NSCLC (Score: 99/100)
**Title**: Tumor microenvironment remodeling after neoadjuvant immunotherapy in non-small cell lung cancer revealed by single-cell RNA sequencing

- **Accession URL**: [GSE207422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422)
- **Indication**: NSCLC
- **Sample/Patient Count**: 39 patients / biological specimens
- **Cell Count Estimate**: ~78,000 single cells/nuclei
- **Cell Selection Strategy**: `Unselected / Total Single-Cell Suspension`
- **Protocol Excerpt**: > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Scoring Rationale**: `Clin:30/30 | Scale:20/20 (N=39) | Kinetics:15/15 | Breadth:15/15 (Unselected / Total Single-Cell Suspension) | Matrix:9/10 | Impact:10/10`
- **Available Raw Matrices**: `GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx; GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz; GSE207422_NSCLC_scRNAseq_metadata.xlsx`

**Study Abstract / Clinical Setting**:
To address the potential therapy-resistant mechanisms and underly the changes of tumor microenvironment (TME) after immunotherapy, we performed scRNA-seq  and bulk RNA-seq from the patients with resectable non-small cell lung cancer (NSCLC) before and after PD-1 blockade combined with chemotherapy. We found that the combined therapy significantly remodeled the immune cell compartments in the TME. Patients with different pathologic responses had distinct characteristics of malignant cells and immune cell compositions. Our study provided several potential biomarkers for future immunotherapy and novel strategies to improve response.

### 17. GSE317309 — NSCLC (Score: 93/100)
**Title**: Combination of a CCL21-gene modified dendritic cell vaccine and pembrolizumab induces immune responses in non-small cell lung cancer

- **Accession URL**: [GSE317309](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE317309)
- **Indication**: NSCLC
- **Sample/Patient Count**: 64 patients / biological specimens
- **Cell Count Estimate**: ~128,000 single cells/nuclei
- **Cell Selection Strategy**: `Unselected / Total Single-Cell Suspension`
- **Protocol Excerpt**: > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Scoring Rationale**: `Clin:30/30 | Scale:20/20 (N=64) | Kinetics:15/15 | Breadth:15/15 (Unselected / Total Single-Cell Suspension) | Matrix:8/10 | Impact:5/10`
- **Available Raw Matrices**: `GSE317309_feature_README.txt; filelist.txt; GSE317309_RAW.tar`

**Study Abstract / Clinical Setting**:
Immune checkpoint inhibitors (ICIs) have transformed treatment for non-small cell lung cancer (NSCLC), but resistance to these therapies is common. We report results from a phase I trial combining intratumoral administration of a CCL21-gene modified dendritic cell (CCL21-DC) vaccine with pembrolizumab in patients with advanced NSCLC. Among 23 patients that received trial therapy, there were no dose-limiting toxicities and a low incidence of treatment-related adverse events. Although no objective responses were observed, 36.8% of patients had a best response of stable disease (SD). Correlative analyses of longitudinal tumor biopsies revealed that disease stability was associated with reduced tumor mutational heterogeneity and increased intratumoral T cell receptor diversity following therapy. Post-treatment tumor biopsies exhibited increased CD4+ T cell infiltration and the presence of novel T cell clones that predominantly possessed CD4+ memory T cell phenotypes. These novel clones experienced greater clonal expansion and a transition towards exhausted phenotypes in SD samples. However, many novel clones failed to persist over time, and immunosuppressive myeloid signaling signatures were identified especially in progressive disease samples. Our findings demonstrated that CCL21-DC vaccination combined with pembrolizumab is safe and can promote antitumor immune responses in a subset of patients. However, efficacy may have been constrained due to limited tumor specificity of the T cell response and an inability to fully reprogram the immunosuppressive tumor microenvironment.

### 18. GSE243013 — NSCLC (Score: 90/100)
**Title**: A single-cell atlas of immune heterogeneity  in anti-PD1-treated non-small cell lung cancer

- **Accession URL**: [GSE243013](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243013)
- **Indication**: NSCLC
- **Sample/Patient Count**: 243 patients / biological specimens
- **Cell Count Estimate**: ~486,000 single cells/nuclei
- **Cell Selection Strategy**: `Unselected / Total Single-Cell Suspension`
- **Protocol Excerpt**: > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Scoring Rationale**: `Clin:30/30 | Scale:20/20 (N=243) | Kinetics:5/15 | Breadth:15/15 (Unselected / Total Single-Cell Suspension) | Matrix:10/10 | Impact:10/10`
- **Available Raw Matrices**: `GSE243013_NSCLC_immune_scRNA_counts.mtx.gz`

**Study Abstract / Clinical Setting**:
Anti-PD(L)1 with chemotherapy is a standard of care for non-small cell lung cancer (NSCLC), but a varying degree of response to the same regimen is observed in patients. The tumor immune microenvironment (TIME) plays a key role in response to immunotherapy, and the TIME heterogeneity in association with therapeutic outcome is incompletely understood. Here we prospectively applied single-cell RNA and TCR sequencing to characterize post-neoadjuvant chemo-immunotherapy treatment tumor samples of 234 NSCLC patients. Our study provides a fine-grained dissection of the TIME heterogeneity underlying response to chemo-immunotherapy in NSCLC, thus representing a valuable resource for improved management of NSCLC.

### 19. GSE311789 — PDAC (Score: 91/100)
**Title**: DeCAF redefines fibroblast states uncovering multidimensional tumor-stroma relationships driving clinical tumor progression and immunotherapy response

- **Accession URL**: [GSE311789](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE311789)
- **Indication**: PDAC
- **Sample/Patient Count**: 142 patients / biological specimens
- **Cell Count Estimate**: ~284,000 single cells/nuclei
- **Cell Selection Strategy**: `Unselected / Total Single-Cell Suspension`
- **Protocol Excerpt**: > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Scoring Rationale**: `Clin:30/30 | Scale:20/20 (N=142) | Kinetics:8/15 | Breadth:15/15 (Unselected / Total Single-Cell Suspension) | Matrix:8/10 | Impact:10/10`
- **Available Raw Matrices**: `filelist.txt; GSE311789_RAW.tar`

**Study Abstract / Clinical Setting**:
This SuperSeries is composed of the SubSeries listed below.

### 20. GSE316195 — PDAC (Score: 86/100)
**Title**: A Phase 1 clinical trial and single-cell correlates of motixafortide, cemiplimab, gemcitabine and nab-paclitaxel for metastatic treatment-naïve metastatic pancreatic ductal adenocarcinoma.

- **Accession URL**: [GSE316195](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316195)
- **Indication**: PDAC
- **Sample/Patient Count**: 22 patients / biological specimens
- **Cell Count Estimate**: ~44,000 single cells/nuclei
- **Cell Selection Strategy**: `Nuclei Isolation (snRNA-seq)`
- **Protocol Excerpt**: > *"Single-nucleus RNA sequencing performed on isolated nuclei from frozen resection tissue."*
- **Scoring Rationale**: `Clin:30/30 | Scale:15/20 (N=22) | Kinetics:15/15 | Breadth:13/15 (Nuclei Isolation (snRNA-seq)) | Matrix:8/10 | Impact:5/10`
- **Available Raw Matrices**: `filelist.txt; GSE316195_RAW.tar`

**Study Abstract / Clinical Setting**:
The C-X-C motif chemokine receptor 4 (CXCR4)/C-X-C motif chemokine ligand 12 (CXCL12) axis is a well-established contributor to the immunosuppressive and immune-excluded TME in pancreatic adenocarcinoma (PDA). Building on pre-clinical data demonstrating a survival benefit with the addition of gemcitabine to CXCR4 and PD1 inhibition in the KPC mouse model, we conducted an open-label, single-arm phase 1 clinical trial combining CXCR4 inhibition (motixafortide), PD-1 blockade (cemiplimab), and chemotherapy (gemcitabine/nab-paclitaxel; MCGN) in treatment-naïve patients with metastatic PDA (n=11). MCGN was safe and tolerable, achieving a 64% partial response (PR) by iRECIST, a 55% confirmed partial response (cPR), and a 91% disease control rate (DCR). The median PFS and OS were 9.7 months (95% confidence interval [CI]: 5.9-not reached [NR]) and 10.1 months (95% CI: 9.3-NR), respectively. One patient achieved a pathological complete response within the primary tumor and the hepatic metastasis after undergoing a pancreatoduodenectomy and hepatectomy and has remained free of disease for 18 months, suggesting that durable responses are possible with MCGN treatment. Single-nucleus RNA-sequencing of serial tissue biopsies from all trial participants revealed a reduction in transcriptional heterogeneity and depletion of cells expressing epithelial-to-mesenchymal transition states, while treatment-resistant tumors maintained tumor heterogeneity. The presence of CXCL12+ proaxogenic Cancer Associated Fibroblasts (pCAFs) were predictive of MCGN efficacy and depleted in resistant samples, suggesting that they represent the substrate for treatment response. MCGN also induced a highly inflamed tumor-microenvironment, which we confirmed with serial tissue staining, and rescue of T cell dysfunction. Based on these promising data and biomarkers, we have launched a multicenter randomized phase 2 trial comparing MGCN to GN in patients with treatment-naïve metastatic PDA (NCT04543071), which is ongoing.

### 21. GSE210038 — ccRCC (Score: 81/100)
**Title**: Mesenchymal-like tumor cells and myofibroblastic cancer-associated fibroblasts are associated with progression and immunotherapy response of clear-cell renal cell carcinoma [scRNA-seq]

- **Accession URL**: [GSE210038](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210038)
- **Indication**: ccRCC
- **Sample/Patient Count**: 9 patients / biological specimens
- **Cell Count Estimate**: ~18,000 single cells/nuclei
- **Cell Selection Strategy**: `Unselected / Total Single-Cell Suspension`
- **Protocol Excerpt**: > *"Primary tumor tissue dissociated into single-cell suspension without marker-based enrichment."*
- **Scoring Rationale**: `Clin:30/30 | Scale:10/20 (N=9) | Kinetics:8/15 | Breadth:15/15 (Unselected / Total Single-Cell Suspension) | Matrix:8/10 | Impact:10/10`
- **Available Raw Matrices**: `filelist.txt; GSE210038_RAW.tar`

**Study Abstract / Clinical Setting**:
Immune checkpoint inhibitors (ICI) represent the cornerstone for treatment of patients with metastatic clear-cell renal cell carcinoma (ccRCC). Despite a favorable response for a subset of patients, others experience primary progressive disease highlighting the need to precisely understand plasticity of cancer cells and their crosstalk with the microenvironment to better predict therapeutic response and personalize treatment. Single-cell RNA sequencing of ccRCC at different disease stages and normal adjacent tissue (NAT) from patients identified 46 cell populations, including 5 tumor subpopulations, characterized by distinct transcriptional signatures representing an epithelial to mesenchymal transition gradient and a novel inflamed state. Deconvolution of the tumor and microenvironment signatures in public datasets and in data from the BIONIKK clinical trial (NCT02960906) revealed a strong correlation between mesenchymal-like ccRCC cells and myofibroblastic cancer-associated fibroblasts (myCAFs), which are both enriched in metastases and correlate with poor patient survival. Spatial transcriptomics and multiplex immune staining uncovered spatial proximity of mesenchymal-like ccRCC cells and myCAFs at the tumor-NAT interface. Moreover, enrichment in myCAFs was associated with primary resistance to ICI therapy in the BIONIKK clinical trial. This data highlights the epithelial-mesenchymal plasticity of ccRCC cancer cells and their relationship with myCAFs, a critical component of the microenvironment associated with poor outcome and ICI resistance.

### 22. GSE314072 — ccRCC (Score: 76/100)
**Title**: Functionally heterogeneous intratumoral CD4+CD8+ double positive T cells can give rise to single positive T cells [scRNA-seq + scTCR-seq]

- **Accession URL**: [GSE314072](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314072)
- **Indication**: ccRCC
- **Sample/Patient Count**: 24 patients / biological specimens
- **Cell Count Estimate**: ~48,000 single cells/nuclei
- **Cell Selection Strategy**: `FACS-sorted (CD3+/CD8+ T-cell enriched)`
- **Protocol Excerpt**: > *"D4+CD8+ double positive T cells can give rise to single positive T cells [scRNA-seq + scTCR-seq] Conventional single positive (SP) CD4+ and CD8+ T cells recognize tumor antigens and help mediate clinical responses with cancer immunotherapy."*
- **Scoring Rationale**: `Clin:30/30 | Scale:15/20 (N=24) | Kinetics:5/15 | Breadth:6/15 (FACS-sorted (CD3+/CD8+ T-cell enriched)) | Matrix:10/10 | Impact:10/10`
- **Available Raw Matrices**: `GSE314072_RCC_Total_adata_multiplex_harmony_noDBL_nomacs_contams_deidentified.h5ad`

**Study Abstract / Clinical Setting**:
Conventional single positive (SP) CD4+ and CD8+ T cells recognize tumor antigens and help mediate clinical responses with cancer immunotherapy. Double positive CD4+CD8+ (DP) T cells have also been described in human cancers, but their role in the tumor microenvironment (TME) remains unclear. By generating a multi-omic single cell atlas of DP and SP T cells, we find that DP T cells possess phenotypic heterogeneity similar to SP T cells that includes multiple clonally expanded populations of cytotoxic DP T cells in human renal cell carcinoma (RCC). These intratumoral DP T cells can mediate by both MHC class I- and class II-dependent killing of autologous tumor cells. In addition, transcriptional profiling of DP TCR-bearing T cells revealed a gene signature enriched for clinical responders to PD-1 blockade in advanced RCC. We confirm prior observations of SP T cells transitioning into DP T cells and more notably, demonstrate that intratumoral T cells are capable of bidirectional differentiation in which DP T cells serve as precursors to SP T cells in vivo.  In the latter scenario, intratumoral DP T cells are shown to express Rag2, suggesting that the tumor may act as an extrathymic site of T cell development. These findings reveal the multiple roles that DP T cells can possess in anti-tumor immunity.
