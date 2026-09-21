# Single-Cell & Spatial Transcriptomics Studies of Immune Checkpoint Therapy Patients

## Reference paper (your starting point)

**Sade-Feldman et al. (2018), *Cell*** — "Defining T cell states associated with response to checkpoint immunotherapy in melanoma." scRNA-seq (Smart-seq2) of 16,291 CD45+ immune cells from 48 tumor samples of 32 metastatic melanoma patients treated with anti-PD-1, anti-CTLA-4 + anti-PD-1, or anti-CTLA-4. Identified two CD8+ T cell states (memory/TCF7+ "CD8_G" vs exhausted "CD8_B") linked to response. Data: GEO **GSE120575**. DOI: [10.1016/j.cell.2018.10.038](https://doi.org/10.1016/j.cell.2018.10.038).

---

## A. Single-cell RNA-seq (scRNA-seq / snRNA-seq / CITE-seq) studies

### 1. Jerby-Arnon et al. (2018), *Cell*
- **Cancer:** Melanoma
- **Therapy:** Anti-PD-1, anti-CTLA-4 (ICB-resistant samples)
- **Modality:** scRNA-seq (Smart-seq2)
- **Patients/samples:** 31 patients, 33 tumors (15 ICB-resistant)
- **Key finding:** Identified a melanoma-intrinsic cell-stress program associated with T cell exclusion and immunotherapy resistance.
- **Data:** GEO **GSE115978** | DOI: [10.1016/j.cell.2018.09.006](https://doi.org/10.1016/j.cell.2018.09.006)

### 2. Li et al. (2019), *Cell*
- **Cancer:** Melanoma
- **Therapy:** Immune checkpoint inhibitors (ICI)
- **Modality:** scRNA-seq (Smart-seq2)
- **Patients/samples:** 29 tumors (untreated or ICI-treated)
- **Key finding:** Characterized immune cell populations and their association with ICI response.
- **Data:** GEO **GSE123139** | DOI: [10.1016/j.cell.2019.08.004](https://doi.org/10.1016/j.cell.2019.08.004)

### 3. de Andrade et al. (2019), *JCI Insight*
- **Cancer:** Melanoma
- **Therapy:** Immune checkpoint inhibitors
- **Modality:** scRNA-seq
- **Patients/samples:** 5 tumors treated with ICI
- **Key finding:** Described immune cell phenotypes in ICI-treated melanoma tumors.
- **Data:** GEO **GSE139249** | DOI: [10.1172/jci.insight.128943](https://doi.org/10.1172/jci.insight.128943)

### 4. Pozniak et al. (2024), *Cell*
- **Cancer:** Melanoma
- **Therapy:** ICB (responders vs non-responders)
- **Modality:** scRNA-seq + spatial multi-omics (10X)
- **Patients/samples:** 28 samples (17 pre-treatment, 11 post-treatment)
- **Key finding:** TCF4-dependent gene regulatory network confers resistance to immunotherapy; targeting TCF4 increases MES-cell immunogenicity and sensitivity to ICB.
- **Data:** KU Leuven RDR, DOI [10.48804/GSAXBN](https://doi.org/10.48804/GSAXBN) | Article: [10.1016/j.cell.2023.11.037](https://doi.org/10.1016/j.cell.2023.11.037)

### 5. Alvarez-Breckenridge et al. (2022), *Cancer Immunology Research*
- **Cancer:** Melanoma brain metastases
- **Therapy:** Immune checkpoint inhibitors
- **Modality:** scRNA-seq (Smart-seq2) + scTCR-seq
- **Patients/samples:** 24 malignant / 25 immune samples (pre/post ICI)
- **Key finding:** Defined microenvironmental correlates of ICI response in melanoma brain metastases.
- **Data:** Single Cell Portal **SCP1493** | DOI: [10.1158/2326-6066.CIR-21-0870](https://doi.org/10.1158/2326-6066.CIR-21-0870)

### 6. Yost et al. (2019), *Nature Medicine*
- **Cancer:** Basal cell carcinoma and squamous cell carcinoma
- **Therapy:** Anti-PD-1
- **Modality:** Paired scRNA-seq + scTCR-seq (10X)
- **Patients/samples:** 79,046 cells from site-matched tumors pre- and post-anti-PD-1
- **Key finding:** PD-1 blockade drives clonal replacement — expanded T cell clones are novel, not reinvigorated pre-existing TILs.
- **Data:** GEO **GSE123813** / **GSE123814** | DOI: [10.1038/s41591-019-0522-3](https://doi.org/10.1038/s41591-019-0522-3)

### 7. Bassez et al. (2021), *Nature Medicine*
- **Cancer:** Breast cancer (TNBC, HER2+, ER+)
- **Therapy:** Anti-PD-1 (pembrolizumab); some patients received neoadjuvant chemo before anti-PD-1
- **Modality:** scRNA-seq + scTCR-seq + CITE-seq (10X)
- **Patients/samples:** 29 anti-PD-1 naive + 11 neoadjuvant chemo then anti-PD-1 (paired pre/post biopsies)
- **Key finding:** Anti-PD-1 induces clonal expansion of PD1+ T cells in ~1/3 of tumors; identified pretreatment immunophenotypes correlated with expansion.
- **Data:** EGA **EGAS00001004809** / **EGAD00001006608**; processed counts at [lambrechtslab.sites.vib.be/en/single-cell](https://lambrechtslab.sites.vib.be/en/single-cell) | DOI: [10.1038/s41591-021-01323-8](https://doi.org/10.1038/s41591-021-01323-8)

### 8. Zhang et al. (2021), *Cancer Cell*
- **Cancer:** Triple-negative breast cancer (TNBC)
- **Therapy:** Atezolizumab (anti-PD-L1) + paclitaxel vs paclitaxel alone
- **Modality:** scRNA-seq (10X)
- **Patients/samples:** 22 patients (half atezolizumab + paclitaxel, half paclitaxel)
- **Key finding:** Identified key immune cell subsets (PD1+ CXCL13+ T cells, IgG+ plasma cells, SPP1+ macrophages) associated with response to PD-L1 blockade.
- **Data:** GEO **GSE169246** | DOI: [10.1016/j.ccell.2021.09.010](https://doi.org/10.1016/j.ccell.2021.09.010)

### 9. Bi et al. (2021), *Cancer Cell*
- **Cancer:** Clear cell renal cell carcinoma (ccRCC)
- **Therapy:** Immune checkpoint blockade (ICI)
- **Modality:** scRNA-seq (10X)
- **Patients/samples:** 6 malignant / 7 immune samples (pre/post ICI; partial responders + stable disease)
- **Key finding:** ICB remodels the RCC microenvironment; characterized cancer and immune cell reprogramming during therapy.
- **Data:** Single Cell Portal **SCP1288**; 3CA: [weizmann.ac.il/sites/3CA/kidney](https://www.weizmann.ac.il/sites/3CA/kidney) | DOI: [10.1016/j.ccell.2021.02.015](https://doi.org/10.1016/j.ccell.2021.02.015)

### 10. Ma et al. (2019), *Cancer Cell*
- **Cancer:** Hepatocellular carcinoma (HCC) and intrahepatic cholangiocarcinoma (iCCA)
- **Therapy:** ICB / immunotherapy
- **Modality:** scRNA-seq (10X)
- **Patients/samples:** HCC: 6 malignant / 9 immune; iCCA: 8 malignant / 10 immune (pre/post treatment)
- **Key finding:** Characterized tumor and immune cell states in liver cancers during immunotherapy.
- **Data:** GEO **GSE125449** | DOI: [10.1016/j.ccell.2019.08.007](https://doi.org/10.1016/j.ccell.2019.08.007)

### 11. Liu et al. (2022), *Nature Cancer*
- **Cancer:** Non-small cell lung cancer (NSCLC)
- **Therapy:** Anti-PD-1 + chemotherapy (combination therapy)
- **Modality:** scRNA-seq + scTCR-seq (10X 5')
- **Patients/samples:** 47 tumor biopsies from 36 patients (33 treatment-naive, 9 post-treatment responsive, 5 post-treatment nonresponsive)
- **Key finding:** Identified precursor exhausted T (Texp) cells that expand in responsive tumors following anti-PD-1; clonal revival rather than reinvigoration of terminal Tex cells.
- **Data:** GEO **GSE179994** | Interactive portal: [nsclcpd1.cancer-pku.cn](http://nsclcpd1.cancer-pku.cn/) | DOI: [10.1038/s43018-021-00292-8](https://doi.org/10.1038/s43018-021-00292-8)

### 12. Liu et al. (2025), *Cell*
- **Cancer:** NSCLC (LUAD + LUSC)
- **Therapy:** Neoadjuvant anti-PD-1 + chemotherapy
- **Modality:** scRNA-seq + scTCR-seq
- **Patients/samples:** 234 NSCLC patients (post-neoadjuvant chemo-immunotherapy surgical samples)
- **Key finding:** Largest single-cell atlas of anti-PD-1-treated NSCLC; fine-grained dissection of TIME heterogeneity and response to chemo-immunotherapy.
- **Data:** GEO **GSE243013** | DOI: [10.1016/j.cell.2025.03.018](https://doi.org/10.1016/j.cell.2025.03.018)

### 13. Luoma et al. (2022), *Cell*
- **Cancer:** Head and neck squamous cell carcinoma (HNSCC / OSCC)
- **Therapy:** Neoadjuvant anti-PD-1 or anti-PD-1 + anti-CTLA-4 (clinical trial NCT02919683)
- **Modality:** scRNA-seq + scTCR-seq
- **Patients/samples:** 29 patients (pre- and post-treatment tumor and blood samples)
- **Key finding:** Tissue-resident memory T cells are early responders to neoadjuvant ICB; PD-1+ KLRG1- CD8+ T cells in pretreatment blood predict pathologic tumor regression.
- **Data:** GEO **GSE200996** | DOI: [10.1016/j.cell.2022.06.018](https://doi.org/10.1016/j.cell.2022.06.018)

### 14. Li et al. (2023), *Cancer Cell*
- **Cancer:** Colorectal cancer (dMMR/MSI-H)
- **Therapy:** Neoadjuvant anti-PD-1 blockade
- **Modality:** scRNA-seq
- **Patients/samples:** 19 dMMR/MSI-H CRC patients receiving neoadjuvant PD-1 blockade
- **Key finding:** In tumors with pathological complete response, concerted decrease in CD8+ Trm-mitotic, CD4+ Tregs, IL1B+ monocytes, and CCL2+ fibroblasts after treatment; identified proinflammatory TME features mediating residual tumor persistence.
- **Data:** GEO **GSE205506** | DOI: [10.1016/j.ccell.2023.04.011](https://doi.org/10.1016/j.ccell.2023.04.011)

### 15. Wu et al. (2023), *BMC Medicine*
- **Cancer:** Metastatic colorectal cancer (MSI-H/dMMR)
- **Therapy:** First-line anti-PD-1 monotherapy
- **Modality:** scRNA-seq (10X)
- **Patients/samples:** 16 MSI-H/dMMR mCRC patients (treatment-resistant vs treatment-sensitive)
- **Key finding:** CD8+ T cells and IL-1β most correlated with anti-PD-1 resistance; IL-1β-driven MDSC infiltration suppresses CD8+ T cells and enhances resistance.
- **Data:** NCBI BioProject **PRJNA932556** (SRA) + Supplementary Table S1 | DOI: [10.1186/s12916-023-02866-y](https://doi.org/10.1186/s12916-023-02866-y)

### 16. Hwang et al. (2022), *Nature Genetics*
- **Cancer:** Pancreatic ductal adenocarcinoma (PDAC)
- **Therapy:** Neoadjuvant FOLFIRINOX ± radiotherapy; 7 patients received losartan and/or nivolumab
- **Modality:** snRNA-seq + spatial transcriptomics (NanoString GeoMx DSP)
- **Patients/samples:** 43 tumors (18 untreated, 25 treated); 224,988 nuclei
- **Key finding:** Defined multicellular communities (classical, squamoid-basaloid, treatment-enriched) and a neural-like progenitor malignant cell program enriched after therapy.
- **Data:** GEO **GSE202051** (snRNA-seq), **GSE199102** (GeoMx); SCP **SCP1089** (untreated), **SCP1096** (treated) | DOI: [10.1038/s41588-022-01134-8](https://doi.org/10.1038/s41588-022-01134-8)

---

## B. Spatial transcriptomics studies

### 17. Meylan et al. (2022), *Immunity*
- **Cancer:** Clear cell renal cell carcinoma (ccRCC)
- **Therapy:** Nivolumab, nivolumab + ipilimumab, or VEGFR-TKI
- **Modality:** Spatial transcriptomics (10X Visium)
- **Patients/samples:** 24 primary tumor sections (12 frozen, 12 FFPE)
- **Key finding:** TLS generate and propagate anti-tumor antibody-producing plasma cells in ccRCC; IgG+ PCs disseminate along fibroblastic tracks into tumor beds.
- **Data:** GEO **GSE175540** | DOI: [10.1016/j.immuni.2022.02.001](https://doi.org/10.1016/j.immuni.2022.02.001)

### 18. Zhang et al. (2024), *Nature Communications*
- **Cancer:** Colorectal cancer (dMMR and pMMR)
- **Therapy:** Neoadjuvant anti-PD-1 antibody (dMMR patients; responders CR/PR and non-responders SD)
- **Modality:** Spatial transcriptomics (Stereo-seq) + scRNA-seq + multiplexed imaging
- **Patients/samples:** 23 CRC patients (16 Stereo-seq samples from 15 patients; 10 scRNA-seq samples from 10 patients); includes 11 anti-PD-1-treated dMMR patients
- **Key finding:** A 300 µm tumor-stroma boundary region regulates immune cell influx; LAMP3+ DCs and CXCL13+ T cells accumulate at the boundary in ICB responders, while CXCL14+ CAFs form a structural barrier in non-responders.
- **Data:** CNCB GSA-Human **PRJCA020107** (controlled access) & [stomics.tech](https://www.stomics.tech/sap/home.html) | DOI: [10.1038/s41467-024-54710-3](https://doi.org/10.1038/s41467-024-54710-3)

### 19. Mebane et al. (2025), *iScience*
- **Cancer:** Triple-negative breast cancer (TNBC)
- **Therapy:** Neoadjuvant pembrolizumab ± SBRT (clinical trial NCT03366844)
- **Modality:** Spatial transcriptomics (NanoString CosMx, ~1,000 genes) + scRNA-seq
- **Patients/samples:** 4 patients (651,683 single cells analyzed)
- **Key finding:** Characterized TLS-like structures and their association with response to pembrolizumab ± radiation.
- **Data:** Zenodo [10.5281/zenodo.14963458](https://doi.org/10.5281/zenodo.14963458); scRNA-seq: GEO **GSE246613** | DOI: [10.1016/j.isci.2025.112808](https://doi.org/10.1016/j.isci.2025.112808)

### 20. Italiano et al. (2022), *Nature Medicine*
- **Cancer:** Soft-tissue sarcoma
- **Therapy:** Pembrolizumab + low-dose cyclophosphamide (PEMBROSARC trial)
- **Modality:** Spatial transcriptomics (NanoString GeoMx DSP — bulk ROI profiling)
- **Patients/samples:** 6 tumors (3 responders, 3 progressive disease)
- **Key finding:** Identified spatial immune features associated with response vs progression on pembrolizumab in sarcoma.
- **Data:** Deidentified data upon request to corresponding author | DOI: [10.1038/s41591-022-01821-3](https://doi.org/10.1038/s41591-022-01821-3)

### 21. Larroquette et al. (2022), *Journal for ImmunoTherapy of Cancer*
- **Cancer:** NSCLC
- **Therapy:** Anti-PD-1/PD-L1
- **Modality:** Spatial transcriptomics (NanoString GeoMx DSP — bulk ROI profiling)
- **Patients/samples:** 152 patients; 16 samples (8 low macrophage infiltration, 8 high)
- **Key finding:** Spatial dissection of macrophage infiltration and its association with ICI response in NSCLC.
- **Data:** Available upon reasonable request to corresponding author | DOI: [10.1136/jitc-2021-003890](https://doi.org/10.1136/jitc-2021-003890)

### 22. Park et al. (2023), *Cancer Research*
- **Cancer:** Gastric cancer
- **Therapy:** ICI
- **Modality:** Spatial transcriptomics (NanoString GeoMx DSP — bulk ROI profiling)
- **Patients/samples:** 12 tumors (5 responders, 7 non-responders)
- **Key finding:** Identified spatial immune features predictive of ICI response in gastric cancer.
- **Data:** AACR 2023 Annual Meeting Abstract (Cancer Res 83(7_Suppl):2262); no public sequencing accession deposited.

### 23. Peyraud et al. (2023 / 2025), *Cancer Research* / *Cell Reports Medicine*
- **Cancer:** NSCLC
- **Therapy:** ICI
- **Modality:** Spatial transcriptomics (NanoString GeoMx DSP — bulk ROI profiling)
- **Patients/samples:** 6 tumors (3 responders, 3 progressive disease)
- **Key finding:** Spatial profiling of immune niches and cell-cell interactions associated with ICI response in NSCLC.
- **Data:** Full paper: *Cell Reports Med* 2025 (6(2):101934, DOI: [10.1016/j.xcrm.2025.101934](https://doi.org/10.1016/j.xcrm.2025.101934)); Data requires French ethics committee approval (CPP du Sud-Ouest et d’Outre-Mer III).

### 24. Liu et al. (2023), *Journal of Hepatology*
- **Cancer:** Hepatocellular carcinoma (HCC)
- **Therapy:** Anti-PD-1
- **Modality:** Spatial transcriptomics (10X Visium spots)
- **Patients/samples:** 8 patients (5 non-responders, 3 responders) + 3 adjacent normal
- **Key finding:** Spatially resolved architecture of the HCC TME associated with anti-PD-1 response.
- **Data:** CNCB GSA-Human (controlled access) | DOI: [10.1016/j.jhep.2023.01.011](https://doi.org/10.1016/j.jhep.2023.01.011)

---

## C. Notable multi-modal / emerging studies

### 25. Krishna et al. (2021), *Cancer Cell*
- **Cancer:** Clear cell renal cell carcinoma (ccRCC)
- **Therapy:** ICI (nivolumab ± ipilimumab)
- **Modality:** scRNA-seq + spatial profiling (multi-regional)
- **Patients/samples:** Multi-regional sampling of ccRCC tumors
- **Key finding:** Linked multiregional immune landscapes and tissue-resident T cells to tumor topology and therapy efficacy.
- **Data:** NCBI BioProject **PRJNA705464** (GEO GSE171306) & [CZ CELLxGENE](https://cellxgene.cziscience.com/collections/3f50314f-bdc9-40c6-8e4a-b0901ebfbe4c) | DOI: [10.1016/j.ccell.2021.03.007](https://doi.org/10.1016/j.ccell.2021.03.007)

### 26. Gondal et al. (2025), *Scientific Data* — Integrated resource
- **Scope:** Compiled 8 scRNA-seq datasets from 9 cancer types (melanoma, BCC, breast, ccRCC, HCC, iCCA, melanoma brain metastasis), 223 patients, 90,270 cancer cells, 265,671 other cells
- **Key value:** Single curated, quality-checked, and explorable resource for ICB-treated patient scRNA-seq data across cancer types with harmonized response annotations.
- **Data:** Zenodo [10.5281/zenodo.10407126](https://doi.org/10.5281/zenodo.10407126) & CZ CELLxGENE
- **Interactive:** [scrnaseqicb.shinyapps.io](https://scrnaseqicb.shinyapps.io/09_shinycell_all/)
- **DOI:** [10.1038/s41597-025-04381-6](https://doi.org/10.1038/s41597-025-04381-6)

---

## D. Borderline / excluded studies

These studies are related but did not fully meet the inclusion criteria of having single-cell or spatial transcriptomics data from ICI-treated patient tumor samples:

| Study | Cancer | Reason for exclusion/borderline |
|---|---|---|
| Tirosh et al. (2016), *Science* | Melanoma | scRNA-seq of treatment-naive metastatic melanoma; no ICI-treated patient samples. Foundational dataset reused in ICB studies. DOI: [10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501) |
| Leader et al. (2021), *Cancer Cell* | NSCLC | scRNA-seq + CITE-seq of treatment-naive NSCLC; LCAM score validated in atezolizumab-treated cohort but primary data are treatment-naive. DOI: [10.1016/j.ccell.2021.10.009](https://doi.org/10.1016/j.ccell.2021.10.009) |
| Zhao et al. (2019), *Nature Medicine* | Glioblastoma | Longitudinal profiling of 66 GBM patients treated with anti-PD-1 (nivolumab/pembrolizumab); includes some scRNA-seq data (4,000 cells from one PTEN-mutated tumor) but primarily bulk RNA-seq/WES. Borderline — scRNA-seq is limited. DOI: [10.1038/s41591-019-0349-y](https://doi.org/10.1038/s41591-019-0349-y) |
| Oh et al. (2020), *GSE149652* | Bladder cancer | scRNA-seq of sorted CD8+ T cells from bladder cancer patients including anti-PD-L1 (atezolizumab)-treated samples. Limited patient scale; data in GEO GSE149652. |
| Zhang et al. (2023, *bioRxiv* preprint) | HCC | Spatial transcriptomics (Visium) of cabozantinib + nivolumab-treated HCC; preprint only, limited verification. |
| Saberzadeh-Ardestani et al. (2025), *Clin Cancer Res* | dMMR CRC | Spatial proteomics (GeoMx nCounter), not transcriptomics; 30 anti-PD-1-treated dMMR CRC patients. DOI: [10.1158/1078-0432.CCR-24-0853](https://doi.org/10.1158/1078-0432.CCR-24-0853) |

---

## Summary table: Systematic Multi-Criteria Checkup

The table below summarizes the checkup across four specific criteria:
1. **Public & Unrestricted?**: Accessible via open repository (GEO, SRA, Zenodo, SCP, CELLxGENE) without Data Access Committee (DAC) approval or formal DTA.
2. **ICB Cohort?**: Contains patient tumor samples treated with immune checkpoint blockade.
3. **scRNA-seq?**: Single-cell or single-nucleus RNA sequencing available (versus spatial spot or bulk ROI profiling).
4. **Objective Response?**: Clinical response annotations (RECIST CR/PR vs SD/PD, or neoadjuvant pCR/MPR) available and mapped.

| # | Study (First author, year) | Cancer | Therapy | Public & Unrestricted? | scRNA-seq? | Objective Response? | Actionable Verdict | Data Accession / Repository | DOI |
|:---:|:---|:---|:---|:---:|:---:|:---:|:---:|:---|:---:|
| 1 | Sade-Feldman (2018) | Melanoma | Anti-PD-1, anti-CTLA-4 | **Yes** | **Yes** | **Yes** (RECIST CR/PR vs SD/PD) | **Tier 1 (Benchmark)** | GEO GSE120575 | [10.1016/j.cell.2018.10.038](https://doi.org/10.1016/j.cell.2018.10.038) |
| 2 | Jerby-Arnon (2018) | Melanoma | Anti-PD-1, anti-CTLA-4 | **Yes** | **Yes** | **Yes** (Prior ICB response) | **Tier 1 (Eligible)** | GEO GSE115978 | [10.1016/j.cell.2018.09.006](https://doi.org/10.1016/j.cell.2018.09.006) |
| 3 | Li (2019) | Melanoma | ICI | **Yes** | **Yes** | **Yes** (Response documented) | **Tier 1 (Eligible)** | GEO GSE123139 | [10.1016/j.cell.2019.08.004](https://doi.org/10.1016/j.cell.2019.08.004) |
| 4 | de Andrade (2019) | Melanoma | ICI | **Yes** | **Yes** | **No (Single Class: PD only)** | **Tier 4 (Ineligible: No Responders)** | GEO GSE139249 | [10.1172/jci.insight.128943](https://doi.org/10.1172/jci.insight.128943) |
| 5 | Pozniak (2024) | Melanoma | ICB | **Yes** | **Yes** | **Yes** (Responders vs NR) | **Tier 1 (Eligible)** | KU Leuven RDR 10.48804/GSAXBN | [10.1016/j.cell.2023.11.037](https://doi.org/10.1016/j.cell.2023.11.037) |
| 6 | Alvarez-Breckenridge (2022) | Melanoma brain mets | ICI | **Yes** | **Yes** | **Yes** (Responders vs NR) | **Tier 1 (Eligible)** | Broad SCP1493 | [10.1158/2326-6066.CIR-21-0870](https://doi.org/10.1158/2326-6066.CIR-21-0870) |
| 7 | Yost (2019) | BCC, SCC | Anti-PD-1 | **Yes** | **Yes** | **Yes** (Tumor regression / R vs NR) | **Tier 1 (Eligible)** | GEO GSE123813 / GSE123814 | [10.1038/s41591-019-0522-3](https://doi.org/10.1038/s41591-019-0522-3) |
| 8 | Bassez (2021) | Breast (TNBC, HER2+, ER+) | Anti-PD-1 (pembro) | **Partial (Processed Open)** | **Yes** | **Surrogate (T-cell expansion E vs NE)** | **Tier 2 (Eligible with Caveat)** | Lambrechts Lab / Zenodo (Raw: EGA EGAS00001004809) | [10.1038/s41591-021-01323-8](https://doi.org/10.1038/s41591-021-01323-8) |
| 9 | Zhang (2021) | TNBC | Atezolizumab + chemo | **Yes** | **Yes** | **Yes** (RECIST CR/PR vs SD/PD) | **Tier 1 (Eligible)** | GEO GSE169246 | [10.1016/j.ccell.2021.09.010](https://doi.org/10.1016/j.ccell.2021.09.010) |
| 10 | Bi (2021) | ccRCC | ICI | **Yes** | **Yes** | **Yes** (PR vs SD) | **Tier 1 (Eligible)** | Broad SCP1288 | [10.1016/j.ccell.2021.02.015](https://doi.org/10.1016/j.ccell.2021.02.015) |
| 11 | Ma (2019) | HCC, iCCA | ICI | **Yes** | **Yes** | **Yes** (RECIST response) | **Tier 1 (Eligible)** | GEO GSE125449 | [10.1016/j.ccell.2019.08.007](https://doi.org/10.1016/j.ccell.2019.08.007) |
| 12 | Liu (2022) | NSCLC | Anti-PD-1 + chemo | **Yes** | **Yes** | **Yes** (Responsive vs NR) | **Tier 1 (Eligible)** | GEO GSE179994 | [10.1038/s43018-021-00292-8](https://doi.org/10.1038/s43018-021-00292-8) |
| 13 | Liu (2025) | NSCLC | Anti-PD-1 + chemo | **Yes** | **Yes** | **Yes (Neoadjuvant MPR/pCR)** | **Tier 2 (Eligible: Neoadjuvant)** | GEO GSE243013 | [10.1016/j.cell.2025.03.018](https://doi.org/10.1016/j.cell.2025.03.018) |
| 14 | Luoma (2022) | HNSCC | Anti-PD-1 ± anti-CTLA-4 | **Yes** | **Yes** | **Yes (Neoadjuvant regression)** | **Tier 2 (Eligible: Neoadjuvant)** | GEO GSE200996 | [10.1016/j.cell.2022.06.018](https://doi.org/10.1016/j.cell.2022.06.018) |
| 15 | Li (2023) | CRC (dMMR/MSI-H) | Anti-PD-1 (neoadjuvant) | **Yes** | **Yes** | **Yes (Neoadjuvant pCR)** | **Tier 2 (Eligible: Neoadjuvant)** | GEO GSE205506 | [10.1016/j.ccell.2023.04.011](https://doi.org/10.1016/j.ccell.2023.04.011) |
| 16 | Wu (2023) | CRC (MSI-H) | Anti-PD-1 | **Yes** | **Yes** | **Yes** (Resistant vs Sensitive) | **Tier 1 (Eligible)** | NCBI PRJNA932556 + Supp S1 | [10.1186/s12916-023-02866-y](https://doi.org/10.1186/s12916-023-02866-y) |
| 17 | Hwang (2022) | PDAC | Nivolumab (subset N=7) | **Yes** | **Yes (snRNA)** | **Partial** | **Tier 2/4 (Small ICB N=7)** | GEO GSE202051, SCP1089/1096 | [10.1038/s41588-022-01134-8](https://doi.org/10.1038/s41588-022-01134-8) |
| 18 | Meylan (2022) | ccRCC | Nivolumab ± ipilimumab | **Yes** | **No (Spatial Visium spots ~55µm)** | **Yes** (RECIST response) | **Tier 3 (Spatial Only)** | GEO GSE175540 | [10.1016/j.immuni.2022.02.001](https://doi.org/10.1016/j.immuni.2022.02.001) |
| 19 | Zhang (2024) | CRC (dMMR/pMMR) | Anti-PD-1 (neoadjuvant) | **Restricted** | **Yes** | **Yes** (CR/PR vs SD) | **Tier 4 (Restricted GSA)** | CNCB GSA PRJCA020107 / stomics.tech | [10.1038/s41467-024-54710-3](https://doi.org/10.1038/s41467-024-54710-3) |
| 20 | Mebane (2025) | TNBC | Pembro ± SBRT | **Yes** | **Yes (CosMx + scRNA)** | **Yes** (Clearance / response) | **Tier 1/2 (Eligible: Small N=4)** | Zenodo 10.5281/zenodo.14963458, GSE246613 | [10.1016/j.isci.2025.112808](https://doi.org/10.1016/j.isci.2025.112808) |
| 21 | Italiano (2022) | Sarcoma | Pembro + cyclophosphamide | **No (Upon request)** | **No (Spatial GeoMx bulk ROI)** | **Yes** (Responders vs PD) | **Tier 4 (Ineligible: Request only)** | Authors upon request (PEMBROSARC) | [10.1038/s41591-022-01821-3](https://doi.org/10.1038/s41591-022-01821-3) |
| 22 | Larroquette (2022) | NSCLC | ICI | **No (Upon request)** | **No (Spatial GeoMx bulk ROI)** | **Yes** (Durable clinical benefit) | **Tier 4 (Ineligible: Request only)** | Authors upon request | [10.1136/jitc-2021-003890](https://doi.org/10.1136/jitc-2021-003890) |
| 23 | Park (2023) | Gastric | ICI | **No (Abstract only)** | **No (Spatial GeoMx bulk ROI)** | **Yes** (5 R vs 7 NR) | **Tier 4 (Ineligible: Abstract only)** | None (AACR Abstract 2262) | AACR 2023 / Cancer Res |
| 24 | Peyraud (2023/2025) | NSCLC | ICI | **No (Ethics approval)** | **No (Spatial GeoMx bulk ROI)** | **Yes** (3 R vs 3 PD) | **Tier 4 (Ineligible: Controlled)** | French Ethics Committee approval | [10.1016/j.xcrm.2025.101934](https://doi.org/10.1016/j.xcrm.2025.101934) |
| 25 | Liu (2023) | HCC | Anti-PD-1 | **Restricted** | **No (Spatial Visium spots ~55µm)** | **Yes** (3 R vs 5 NR) | **Tier 3/4 (Spatial Only / GSA)** | CNCB GSA-Human | [10.1016/j.jhep.2023.01.011](https://doi.org/10.1016/j.jhep.2023.01.011) |
| 26 | Krishna (2021) | ccRCC | ICI (nivo ± ipi) | **Yes** | **Yes** | **Yes** (Therapy efficacy) | **Tier 1 (Eligible)** | BioProject PRJNA705464 / CZ CELLxGENE | [10.1016/j.ccell.2021.03.007](https://doi.org/10.1016/j.ccell.2021.03.007) |
| 27 | Gondal (2025) | 9 cancer types (resource) | ICB | **Yes** | **Yes** | **Yes** (Harmonized R vs NR) | **Tier 1 (Meta-Resource)** | Zenodo 10.5281/zenodo.10407126 / CELLxGENE | [10.1038/s41597-025-04381-6](https://doi.org/10.1038/s41597-025-04381-6) |
