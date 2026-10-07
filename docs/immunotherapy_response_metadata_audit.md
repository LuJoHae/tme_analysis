# Immunotherapy Response Metadata Audit Report

This document provides the empirical verification results for all **350 single-cell cohorts** (Tier 0 Benchmark Core, Tier 1 ICB Response / Treated, and on-disk candidates) inspected on server `olm`.

## Executive Summary Statistics
- **Total Cohorts Audited**: 350
- **Verified Response Ground Truth (Class A + B)**: **8 cohorts (2.3%)**
  - **Class A (Direct Cell-Level Response)**: 3 cohorts
  - **Class B (Patient-Level Table / GEO Matrix)**: 5 cohorts
- **Class C (ICB Treated Only, No Response Breakdown)**: 28 cohorts
- **Class D (False Positive / Pre-Clinical / In Vitro)**: 314 cohorts

## Master Audit Table

| Accession | Indication | Original Tier | Veracity Class | Verified? | Response Column | Patient ID Key | Unique Response Categories | Source File |
| :--- | :--- | :--- | :--- | :---: | :--- | :--- | :--- | :--- |
| **GSE344166** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE300446** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE294273** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE320040** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE286410** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE303948** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE270464** | Melanoma | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **GSE256291** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE244983** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE210963** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE242477** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE218429** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE198265** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE338555** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE324655** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE317349** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE230574** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE269936** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE317309** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE276139** | NSCLC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **GSE253718** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE233203** | NSCLC | Tier 1 (ICB Response) | **Class B** | ✓ YES | `therapeutic response` | `!Sample_geo_accession` | `Non-response, Response` | GEO Series Matrix (Sample Characteristics) |
| **GSE270148** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE285888** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE243013** | NSCLC | Tier 0 (Benchmark Core) | **Class B** | ✓ YES | `pathological_response` | `sampleID` | `non-MPR` | GSE243013_NSCLC_immune_scRNA_metadata.csv.gz |
| **GSE241934** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE223779** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE207422** | NSCLC | Tier 0 (Benchmark Core) | **Class B** | ✓ YES | `Pathologic Response` | `Sample` | `MPR, MPR (pCR), NMPR` | GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx (sheet: sheet1) |
| **GSE216069** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE192402** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE211068** | Melanoma | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **GSE346380** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE339453** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE333596** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE327167** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE311609** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE305872** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE316782** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE307811** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE308103** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE303762** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE294109** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE289672** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE314072** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE285701** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE272610** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE220313** | ccRCC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **GSE210038** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE223808** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE328692** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE304262** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE254444** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE254498** | ccRCC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **GSE304466** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE310802** | Bladder | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE302781** | Bladder | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE326225** | Bladder | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE326854** | Bladder | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE301651** | Bladder | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE277524** | Bladder | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE250523** | Bladder | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE267718** | Bladder | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE222315** | Bladder | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE183556** | Bladder | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE176249** | Bladder | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE172433** | Bladder | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE192575** | Bladder | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE302453** | Breast | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **GSE300475** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE229723** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE274141** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE274139** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE262288** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE254991** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE199219** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE246613** | Breast | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **GSE222859** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE212707** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE329389** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE303346** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE337706** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE299267** | Breast | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **GSE325982** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE332708** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE331487** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE309616** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE281490** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE281488** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE300628** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE230327** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE278406** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE236581** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE235917** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE205506** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE216534** | CRC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **GSE164522** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE188711** | CRC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **GSE146771** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE271690** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE294300** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE336564** | CRC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **GSE312260** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE335811** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE311338** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE330797** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE315534** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE312804** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE309346** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE299651** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE270767** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE274321** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE296954** | HNSCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE296867** | HNSCC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **GSE301720** | HNSCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE301741** | HNSCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE287301** | HNSCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE247582** | HNSCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE200996** | HNSCC | Tier 1 (ICB Response) | **Class B** | ✓ YES | `Path_response` | `Patient_ID` | `High, Medium` | GSE200996_CD4.tumor.single.cell.meta.data.txt.gz |
| **GSE339480** | HNSCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE268014** | HNSCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE327189** | HNSCC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **GSE296771** | HNSCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE322620** | HNSCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE310797** | HNSCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE198315** | HNSCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE280982** | HNSCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE286935** | HNSCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE281978** | HNSCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE270680** | Gastric | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE239676** | Gastric | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE234209** | Gastric | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE321676** | Gastric | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE275648** | Gastric | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE308231** | Gastric | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE246662** | Gastric | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE232733** | Gastric | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE228598** | Gastric | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE168537** | Gastric | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE212212** | Gastric | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE184198** | Gastric | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE112302** | Gastric | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE319709** | HCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE318418** | HCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE318420** | HCC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **GSE313642** | HCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE278324** | HCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE281110** | HCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE265770** | HCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE272348** | HCC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **GSE272347** | HCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE255830** | HCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE245906** | HCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE233405** | HCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE224411** | HCC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **GSE215428** | HCC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **GSE291757** | HCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE326201** | HCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE320155** | HCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE208308** | HCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE208307** | HCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE299340** | HCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE282343** | HCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE290925** | HCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE335452** | PDAC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **GSE311789** | PDAC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE311788** | PDAC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE283206** | PDAC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE279781** | PDAC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE205354** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE205049** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE212966** | PDAC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE211644** | PDAC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **GSE202051** | PDAC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE156405** | PDAC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **GSE348275** | PDAC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **GSE348038** | PDAC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE347847** | PDAC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE316195** | PDAC | Tier 0 (Benchmark Core) | **Class B** | ✓ YES | `response` | `!Sample_geo_accession` | `PD, PR, SD` | GEO Series Matrix (Sample Characteristics) |
| **GSE327056** | PDAC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE318413** | PDAC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **GSE300154** | PDAC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE312209** | PDAC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE284392** | PDAC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE288067** | PDAC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE160977** | PDAC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **GSE291124** | PDAC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_05a8c945** | CRC | Tier 1 (ICB Response) | **Class A** | ✓ YES | `RECIST` | `ENA_sample_accession` | `CR: complete response, NE: inevaluable, PD: progressive disease, PR: partial response, SD: stable disease` | 05a8c945-bc12-414f-960d-a31943bbcdd1.h5ad |
| **CELLxGENE_ca140407** | Gastric | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_6f9de485** | Breast | Tier 1 (ICB Response) | **Class A** | ✓ YES | `pCR_status` | `donor_id` | `Excluded, RD, pCR` | 6f9de485-58cd-4342-bfc4-b3d3dd223aa8.h5ad |
| **CELLxGENE_f0e0575d** | Bladder | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_e2094676** | Bladder | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_7ee4b15b** | Bladder | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_670a9f65** | Bladder | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_024581e3** | Bladder | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_714e6bc2** | HNSCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_624d92e2** | HNSCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_60acb72d** | HNSCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_01ff5cf0** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_e3ed2ba4** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_ef7bb7f0** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_e6aaf5a4** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_dc6b1e06** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_9c235282** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_7be23e52** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_763d1d88** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_40a0ade8** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_278eac3f** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_1d54fb17** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_19053a82** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_0d3807bf** | Gastric | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_ef0d813e** | CRC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **CELLxGENE_829a3cd1** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_5ee552f5** | CRC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **CELLxGENE_4b5afdf9** | CRC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **CELLxGENE_387acac5** | CRC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **CELLxGENE_2e95d453** | CRC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **CELLxGENE_2554a654** | CRC | Tier 1 (ICB Treated) | **Class C** | ✗ NO | `` | `` | `None` | Files present, no outcome breakdown |
| **CELLxGENE_ed880090** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_de5416ef** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_75011e96** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_5a9cfb44** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_ee141ea4** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_e5c614b8** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_dd1913a6** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_c829c294** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_b0d9408e** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_9adb1b29** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_5d2c013d** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_540e4c1a** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_480f9371** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_2cc628d1** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_10bb68cf** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_7b20c613** | Melanoma | Tier 1 (ICB Response) | **Class A** | ✓ YES | `Combined_outcome` | `PMID_donor_id` | `Favourable, UT, Unfavourable, n/a` | 7b20c613-9add-43d1-87e9-defd3d9b9f8c.h5ad |
| **CELLxGENE_68b6114f** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_fbdd8c17** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_1e4214ce** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_55ca4411** | HNSCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_ff4cfa86** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_f12ab0e6** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_ec423499** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_e2824739** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_e06e9bf3** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_cd6398a9** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_c7d0def0** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_aa6f371d** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_a6c0143c** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_a6347c54** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_9c5f68fc** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_9237e573** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_7432b873** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_71e44b30** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_71513028** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_6c87755e** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_6384d8b8** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_59d14a35** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_54d56674** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_494faa16** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_48b55b2b** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_44941fdb** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_39f6fec9** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_34f5307e** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_2dd73feb** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_1884e651** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_1637e817** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_12c868c6** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_0f9d1892** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_05a49baa** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_02aa7750** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_f354e4c3** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_aafb780d** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_a6b0f655** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_80466231** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_6f0858c0** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_24dbd26d** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_f7af19e4** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_c0d43178** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_b5753bee** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_a73f7983** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_879bb6df** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_7fe57023** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_7ba1a805** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_74e80fd1** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_729f397a** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_2d821164** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_297b5b89** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_2916b663** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_1e191a00** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_15b98664** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_b617ee1b** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_f25a532c** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_c5ac3ec2** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_c3f74413** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_a6046b15** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_a45125aa** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_81328f3f** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_75548d10** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_60ac2657** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_53d62b10** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_30437616** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_27cd5ac9** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_252438d3** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_24c31c8c** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_104cfa2a** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_0671c0d4** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_05f813a4** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_2f05ab20** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_e500acbf** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_5d3fc988** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_3e4e2c8e** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_9fddb063** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_933497dc** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_7357bdd2** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_4cdd25a4** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_11a3244a** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_0c86f0de** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_04d87de6** | Breast | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_f339fe89** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_dcb7c544** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_d4dc4cfe** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_cc43509e** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_b3052902** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_ae4552dc** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_9b1437bb** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_89972213** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_76bb43ff** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_6d243918** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_50c4a6d6** | Melanoma | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_4aefec71** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_32ffc3a7** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_18fd0190** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_07efa1c3** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_02faf712** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_eaf0c852** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_5af90777** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_318cb2a6** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_bd65a70f** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_9f222629** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_d41f45c1** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_be39785b** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_4c6f9f26** | ccRCC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_7bb64315** | Gastric | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_d6dfdef1** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_6a270451** | CRC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_232f6a5a** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_1e6a6ef9** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_f64e1be1** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_e9175006** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_d4cfefa0** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_d224c8e0** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |
| **CELLxGENE_a6858c10** | NSCLC | Tier 2 (Baseline Atlas) | **Class D** | ✗ NO | `` | `` | `None` | None |

---

## Detailed Cohort Profiles & Notes

### 1. GSE344166 — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 2. GSE300446 — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 3. GSE294273 — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 4. GSE320040 — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 5. GSE286410 — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 6. GSE303948 — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 7. GSE270464 — Melanoma (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 8. GSE256291 — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 9. GSE244983 — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 10. GSE210963 — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 11. GSE242477 — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 12. GSE218429 — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 13. GSE198265 — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 14. GSE338555 — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 15. GSE324655 — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 16. GSE317349 — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 17. GSE230574 — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 18. GSE269936 — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 19. GSE317309 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 20. GSE276139 — NSCLC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 21. GSE253718 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 22. GSE233203 — NSCLC (Class B)
**Verification Status**: `VERIFIED RESPONSE GROUND TRUTH`

- **Original Tier**: Tier 1 (ICB Response)
- **Recommended Tier**: **Tier 1 (ICB Response)**
- **Veracity Class**: `Class B`
- **Endpoint Type**: `RECIST_GEO`
- **Response Column**: `therapeutic response`
- **Patient / Sample Key**: `!Sample_geo_accession`
- **Response Categories Found**: `['Non-response', 'Response']`
- **Metadata Source**: `GEO Series Matrix (Sample Characteristics)`
- **Audit Notes**: Verified in GEO Sample Characteristics (therapeutic response: ('Non-response', 'Response'))

### 23. GSE270148 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 24. GSE285888 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 25. GSE243013 — NSCLC (Class B)
**Verification Status**: `VERIFIED RESPONSE GROUND TRUTH`

- **Original Tier**: Tier 0 (Benchmark Core)
- **Recommended Tier**: **Tier 0 (Benchmark Core)**
- **Veracity Class**: `Class B`
- **Endpoint Type**: `RECIST`
- **Response Column**: `pathological_response`
- **Patient / Sample Key**: `sampleID`
- **Response Categories Found**: `['non-MPR']`
- **Metadata Source**: `GSE243013_NSCLC_immune_scRNA_metadata.csv.gz`
- **Audit Notes**: Verified in metadata table (GSE243013_NSCLC_immune_scRNA_metadata.csv.gz -> col: pathological_response, values: ('non-MPR',))

### 26. GSE241934 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 27. GSE223779 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 28. GSE207422 — NSCLC (Class B)
**Verification Status**: `VERIFIED RESPONSE GROUND TRUTH`

- **Original Tier**: Tier 0 (Benchmark Core)
- **Recommended Tier**: **Tier 0 (Benchmark Core)**
- **Veracity Class**: `Class B`
- **Endpoint Type**: `RECIST`
- **Response Column**: `Pathologic Response`
- **Patient / Sample Key**: `Sample`
- **Response Categories Found**: `['MPR', 'MPR (pCR)', 'NMPR']`
- **Metadata Source**: `GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx (sheet: sheet1)`
- **Audit Notes**: Verified in clinical Excel supplement (GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx (sheet: sheet1) -> col: Pathologic Response, values: ('MPR', 'MPR (pCR)', 'NMPR'))

### 29. GSE216069 — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 30. GSE192402 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 31. GSE211068 — Melanoma (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 32. GSE346380 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 33. GSE339453 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 34. GSE333596 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 35. GSE327167 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 36. GSE311609 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 37. GSE305872 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 38. GSE316782 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 39. GSE307811 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 40. GSE308103 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 41. GSE303762 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 42. GSE294109 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 43. GSE289672 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 44. GSE314072 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 45. GSE285701 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 46. GSE272610 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 47. GSE220313 — ccRCC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 48. GSE210038 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 49. GSE223808 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 50. GSE328692 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 51. GSE304262 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 52. GSE254444 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 53. GSE254498 — ccRCC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 54. GSE304466 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 55. GSE310802 — Bladder (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 56. GSE302781 — Bladder (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 57. GSE326225 — Bladder (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 58. GSE326854 — Bladder (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 59. GSE301651 — Bladder (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 60. GSE277524 — Bladder (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 61. GSE250523 — Bladder (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 62. GSE267718 — Bladder (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 63. GSE222315 — Bladder (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 64. GSE183556 — Bladder (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 65. GSE176249 — Bladder (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 66. GSE172433 — Bladder (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 67. GSE192575 — Bladder (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 68. GSE302453 — Breast (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 69. GSE300475 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 70. GSE229723 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 71. GSE274141 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 72. GSE274139 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 73. GSE262288 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 74. GSE254991 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 75. GSE199219 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 76. GSE246613 — Breast (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 77. GSE222859 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 78. GSE212707 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 79. GSE329389 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 80. GSE303346 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 81. GSE337706 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 82. GSE299267 — Breast (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 83. GSE325982 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 84. GSE332708 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 85. GSE331487 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 86. GSE309616 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 87. GSE281490 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 88. GSE281488 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 89. GSE300628 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 90. GSE230327 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 91. GSE278406 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 92. GSE236581 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 93. GSE235917 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 94. GSE205506 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 95. GSE216534 — CRC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 96. GSE164522 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 97. GSE188711 — CRC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 98. GSE146771 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 99. GSE271690 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 100. GSE294300 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 101. GSE336564 — CRC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 102. GSE312260 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 103. GSE335811 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 104. GSE311338 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 105. GSE330797 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 106. GSE315534 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 107. GSE312804 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 108. GSE309346 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 109. GSE299651 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 110. GSE270767 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 111. GSE274321 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 112. GSE296954 — HNSCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 113. GSE296867 — HNSCC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 114. GSE301720 — HNSCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 115. GSE301741 — HNSCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 116. GSE287301 — HNSCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 117. GSE247582 — HNSCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 118. GSE200996 — HNSCC (Class B)
**Verification Status**: `VERIFIED RESPONSE GROUND TRUTH`

- **Original Tier**: Tier 1 (ICB Response)
- **Recommended Tier**: **Tier 1 (ICB Response)**
- **Veracity Class**: `Class B`
- **Endpoint Type**: `RECIST`
- **Response Column**: `Path_response`
- **Patient / Sample Key**: `Patient_ID`
- **Response Categories Found**: `['High', 'Medium']`
- **Metadata Source**: `GSE200996_CD4.tumor.single.cell.meta.data.txt.gz`
- **Audit Notes**: Verified in metadata table (GSE200996_CD4.tumor.single.cell.meta.data.txt.gz -> col: Path_response, values: ('High', 'Medium'))

### 119. GSE339480 — HNSCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 120. GSE268014 — HNSCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 121. GSE327189 — HNSCC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 122. GSE296771 — HNSCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 123. GSE322620 — HNSCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 124. GSE310797 — HNSCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 125. GSE198315 — HNSCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 126. GSE280982 — HNSCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 127. GSE286935 — HNSCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 128. GSE281978 — HNSCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 129. GSE270680 — Gastric (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 130. GSE239676 — Gastric (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 131. GSE234209 — Gastric (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 132. GSE321676 — Gastric (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 133. GSE275648 — Gastric (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 134. GSE308231 — Gastric (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 135. GSE246662 — Gastric (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 136. GSE232733 — Gastric (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 137. GSE228598 — Gastric (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 138. GSE168537 — Gastric (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 139. GSE212212 — Gastric (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 140. GSE184198 — Gastric (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 141. GSE112302 — Gastric (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 142. GSE319709 — HCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 143. GSE318418 — HCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 144. GSE318420 — HCC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 145. GSE313642 — HCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 146. GSE278324 — HCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 147. GSE281110 — HCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 148. GSE265770 — HCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 149. GSE272348 — HCC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 150. GSE272347 — HCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 151. GSE255830 — HCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 152. GSE245906 — HCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 153. GSE233405 — HCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 154. GSE224411 — HCC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 155. GSE215428 — HCC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 156. GSE291757 — HCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 157. GSE326201 — HCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 158. GSE320155 — HCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 159. GSE208308 — HCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 160. GSE208307 — HCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 161. GSE299340 — HCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 162. GSE282343 — HCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 163. GSE290925 — HCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 164. GSE335452 — PDAC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 165. GSE311789 — PDAC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 166. GSE311788 — PDAC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 167. GSE283206 — PDAC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 168. GSE279781 — PDAC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 169. GSE205354 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 170. GSE205049 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 171. GSE212966 — PDAC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 172. GSE211644 — PDAC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 173. GSE202051 — PDAC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 174. GSE156405 — PDAC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 175. GSE348275 — PDAC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 176. GSE348038 — PDAC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 177. GSE347847 — PDAC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 178. GSE316195 — PDAC (Class B)
**Verification Status**: `VERIFIED RESPONSE GROUND TRUTH`

- **Original Tier**: Tier 0 (Benchmark Core)
- **Recommended Tier**: **Tier 0 (Benchmark Core)**
- **Veracity Class**: `Class B`
- **Endpoint Type**: `RECIST_GEO`
- **Response Column**: `response`
- **Patient / Sample Key**: `!Sample_geo_accession`
- **Response Categories Found**: `['PD', 'PR', 'SD']`
- **Metadata Source**: `GEO Series Matrix (Sample Characteristics)`
- **Audit Notes**: Verified in GEO Sample Characteristics (response: ('PD', 'PR', 'SD'))

### 179. GSE327056 — PDAC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 180. GSE318413 — PDAC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 181. GSE300154 — PDAC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 182. GSE312209 — PDAC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 183. GSE284392 — PDAC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 184. GSE288067 — PDAC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 185. GSE160977 — PDAC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 186. GSE291124 — PDAC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 187. CELLxGENE_05a8c945 — CRC (Class A)
**Verification Status**: `VERIFIED RESPONSE GROUND TRUTH`

- **Original Tier**: Tier 1 (ICB Response)
- **Recommended Tier**: **Tier 1 (ICB Response)**
- **Veracity Class**: `Class A`
- **Endpoint Type**: `RECIST`
- **Response Column**: `RECIST`
- **Patient / Sample Key**: `ENA_sample_accession`
- **Response Categories Found**: `['CR: complete response', 'NE: inevaluable', 'PD: progressive disease', 'PR: partial response', 'SD: stable disease']`
- **Metadata Source**: `05a8c945-bc12-414f-960d-a31943bbcdd1.h5ad`
- **Audit Notes**: Verified directly in cell AnnData .obs table (RECIST: ('CR: complete response', 'NE: inevaluable', 'PD: progressive disease', 'PR: partial response', 'SD: stable disease'))

### 188. CELLxGENE_ca140407 — Gastric (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 189. CELLxGENE_6f9de485 — Breast (Class A)
**Verification Status**: `VERIFIED RESPONSE GROUND TRUTH`

- **Original Tier**: Tier 1 (ICB Response)
- **Recommended Tier**: **Tier 1 (ICB Response)**
- **Veracity Class**: `Class A`
- **Endpoint Type**: `RECIST`
- **Response Column**: `pCR_status`
- **Patient / Sample Key**: `donor_id`
- **Response Categories Found**: `['Excluded', 'RD', 'pCR']`
- **Metadata Source**: `6f9de485-58cd-4342-bfc4-b3d3dd223aa8.h5ad`
- **Audit Notes**: Verified directly in cell AnnData .obs table (pCR_status: ('Excluded', 'RD', 'pCR'))

### 190. CELLxGENE_f0e0575d — Bladder (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 191. CELLxGENE_e2094676 — Bladder (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 192. CELLxGENE_7ee4b15b — Bladder (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 193. CELLxGENE_670a9f65 — Bladder (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 194. CELLxGENE_024581e3 — Bladder (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 195. CELLxGENE_714e6bc2 — HNSCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 196. CELLxGENE_624d92e2 — HNSCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 197. CELLxGENE_60acb72d — HNSCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 198. CELLxGENE_01ff5cf0 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 199. CELLxGENE_e3ed2ba4 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 200. CELLxGENE_ef7bb7f0 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 201. CELLxGENE_e6aaf5a4 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 202. CELLxGENE_dc6b1e06 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 203. CELLxGENE_9c235282 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 204. CELLxGENE_7be23e52 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 205. CELLxGENE_763d1d88 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 206. CELLxGENE_40a0ade8 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 207. CELLxGENE_278eac3f — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 208. CELLxGENE_1d54fb17 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 209. CELLxGENE_19053a82 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 210. CELLxGENE_0d3807bf — Gastric (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 211. CELLxGENE_ef0d813e — CRC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 212. CELLxGENE_829a3cd1 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 213. CELLxGENE_5ee552f5 — CRC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 214. CELLxGENE_4b5afdf9 — CRC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 215. CELLxGENE_387acac5 — CRC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 216. CELLxGENE_2e95d453 — CRC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 217. CELLxGENE_2554a654 — CRC (Class C)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 1 (ICB Treated)
- **Recommended Tier**: **Tier 1 (ICB Treated)**
- **Veracity Class**: `Class C`
- **Endpoint Type**: `Setting_Only`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `Files present, no outcome breakdown`
- **Audit Notes**: Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.

### 218. CELLxGENE_ed880090 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 219. CELLxGENE_de5416ef — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 220. CELLxGENE_75011e96 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 221. CELLxGENE_5a9cfb44 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 222. CELLxGENE_ee141ea4 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 223. CELLxGENE_e5c614b8 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 224. CELLxGENE_dd1913a6 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 225. CELLxGENE_c829c294 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 226. CELLxGENE_b0d9408e — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 227. CELLxGENE_9adb1b29 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 228. CELLxGENE_5d2c013d — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 229. CELLxGENE_540e4c1a — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 230. CELLxGENE_480f9371 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 231. CELLxGENE_2cc628d1 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 232. CELLxGENE_10bb68cf — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 233. CELLxGENE_7b20c613 — Melanoma (Class A)
**Verification Status**: `VERIFIED RESPONSE GROUND TRUTH`

- **Original Tier**: Tier 1 (ICB Response)
- **Recommended Tier**: **Tier 1 (ICB Response)**
- **Veracity Class**: `Class A`
- **Endpoint Type**: `RECIST`
- **Response Column**: `Combined_outcome`
- **Patient / Sample Key**: `PMID_donor_id`
- **Response Categories Found**: `['Favourable', 'UT', 'Unfavourable', 'n/a']`
- **Metadata Source**: `7b20c613-9add-43d1-87e9-defd3d9b9f8c.h5ad`
- **Audit Notes**: Verified directly in cell AnnData .obs table (Combined_outcome: ('Favourable', 'UT', 'Unfavourable', 'n/a'))

### 234. CELLxGENE_68b6114f — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 235. CELLxGENE_fbdd8c17 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 236. CELLxGENE_1e4214ce — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 237. CELLxGENE_55ca4411 — HNSCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 238. CELLxGENE_ff4cfa86 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 239. CELLxGENE_f12ab0e6 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 240. CELLxGENE_ec423499 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 241. CELLxGENE_e2824739 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 242. CELLxGENE_e06e9bf3 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 243. CELLxGENE_cd6398a9 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 244. CELLxGENE_c7d0def0 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 245. CELLxGENE_aa6f371d — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 246. CELLxGENE_a6c0143c — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 247. CELLxGENE_a6347c54 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 248. CELLxGENE_9c5f68fc — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 249. CELLxGENE_9237e573 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 250. CELLxGENE_7432b873 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 251. CELLxGENE_71e44b30 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 252. CELLxGENE_71513028 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 253. CELLxGENE_6c87755e — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 254. CELLxGENE_6384d8b8 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 255. CELLxGENE_59d14a35 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 256. CELLxGENE_54d56674 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 257. CELLxGENE_494faa16 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 258. CELLxGENE_48b55b2b — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 259. CELLxGENE_44941fdb — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 260. CELLxGENE_39f6fec9 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 261. CELLxGENE_34f5307e — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 262. CELLxGENE_2dd73feb — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 263. CELLxGENE_1884e651 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 264. CELLxGENE_1637e817 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 265. CELLxGENE_12c868c6 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 266. CELLxGENE_0f9d1892 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 267. CELLxGENE_05a49baa — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 268. CELLxGENE_02aa7750 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 269. CELLxGENE_f354e4c3 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 270. CELLxGENE_aafb780d — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 271. CELLxGENE_a6b0f655 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 272. CELLxGENE_80466231 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 273. CELLxGENE_6f0858c0 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 274. CELLxGENE_24dbd26d — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 275. CELLxGENE_f7af19e4 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 276. CELLxGENE_c0d43178 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 277. CELLxGENE_b5753bee — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 278. CELLxGENE_a73f7983 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 279. CELLxGENE_879bb6df — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 280. CELLxGENE_7fe57023 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 281. CELLxGENE_7ba1a805 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 282. CELLxGENE_74e80fd1 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 283. CELLxGENE_729f397a — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 284. CELLxGENE_2d821164 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 285. CELLxGENE_297b5b89 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 286. CELLxGENE_2916b663 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 287. CELLxGENE_1e191a00 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 288. CELLxGENE_15b98664 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 289. CELLxGENE_b617ee1b — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 290. CELLxGENE_f25a532c — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 291. CELLxGENE_c5ac3ec2 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 292. CELLxGENE_c3f74413 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 293. CELLxGENE_a6046b15 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 294. CELLxGENE_a45125aa — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 295. CELLxGENE_81328f3f — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 296. CELLxGENE_75548d10 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 297. CELLxGENE_60ac2657 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 298. CELLxGENE_53d62b10 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 299. CELLxGENE_30437616 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 300. CELLxGENE_27cd5ac9 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 301. CELLxGENE_252438d3 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 302. CELLxGENE_24c31c8c — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 303. CELLxGENE_104cfa2a — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 304. CELLxGENE_0671c0d4 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 305. CELLxGENE_05f813a4 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 306. CELLxGENE_2f05ab20 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 307. CELLxGENE_e500acbf — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 308. CELLxGENE_5d3fc988 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 309. CELLxGENE_3e4e2c8e — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 310. CELLxGENE_9fddb063 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 311. CELLxGENE_933497dc — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 312. CELLxGENE_7357bdd2 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 313. CELLxGENE_4cdd25a4 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 314. CELLxGENE_11a3244a — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 315. CELLxGENE_0c86f0de — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 316. CELLxGENE_04d87de6 — Breast (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 317. CELLxGENE_f339fe89 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 318. CELLxGENE_dcb7c544 — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 319. CELLxGENE_d4dc4cfe — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 320. CELLxGENE_cc43509e — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 321. CELLxGENE_b3052902 — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 322. CELLxGENE_ae4552dc — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 323. CELLxGENE_9b1437bb — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 324. CELLxGENE_89972213 — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 325. CELLxGENE_76bb43ff — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 326. CELLxGENE_6d243918 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 327. CELLxGENE_50c4a6d6 — Melanoma (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 328. CELLxGENE_4aefec71 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 329. CELLxGENE_32ffc3a7 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 330. CELLxGENE_18fd0190 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 331. CELLxGENE_07efa1c3 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 332. CELLxGENE_02faf712 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 333. CELLxGENE_eaf0c852 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 334. CELLxGENE_5af90777 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 335. CELLxGENE_318cb2a6 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 336. CELLxGENE_bd65a70f — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 337. CELLxGENE_9f222629 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 338. CELLxGENE_d41f45c1 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 339. CELLxGENE_be39785b — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 340. CELLxGENE_4c6f9f26 — ccRCC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 341. CELLxGENE_7bb64315 — Gastric (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 342. CELLxGENE_d6dfdef1 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 343. CELLxGENE_6a270451 — CRC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 344. CELLxGENE_232f6a5a — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 345. CELLxGENE_1e6a6ef9 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 346. CELLxGENE_f64e1be1 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 347. CELLxGENE_e9175006 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 348. CELLxGENE_d4cfefa0 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 349. CELLxGENE_d224c8e0 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.

### 350. CELLxGENE_a6858c10 — NSCLC (Class D)
**Verification Status**: `UNVERIFIED / SETTING ONLY`

- **Original Tier**: Tier 2 (Baseline Atlas)
- **Recommended Tier**: **Tier 2 (Baseline Atlas)**
- **Veracity Class**: `Class D`
- **Endpoint Type**: `None`
- **Response Column**: ``
- **Patient / Sample Key**: ``
- **Response Categories Found**: `[]`
- **Metadata Source**: `None`
- **Audit Notes**: No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.
