= 1. Introduction

The clinical efficacy of Immune Checkpoint Inhibitors (ICIs) is governed by a complex interplay between cancer-intrinsic alterations and cancer-extrinsic microenvironmental features @addala2024computational. Traditional clinical biomarkers—such as PD-L1 expression by immunohistochemistry or Tumor Mutational Burden (TMB) by exome sequencing—often function as imperfect predictors when evaluated in isolation, owing to variable antibody assays, cohort-specific positive thresholds, and differences in tumor purity @addala2024computational. Extensive meta-analyses have demonstrated that TMB, PD-L1 status, and transcriptional signatures represent largely independent predictors of immunotherapy response, suggesting that multi-layered, integrative models can provide substantial diagnostic improvements @addala2024computational.

To characterize these multi-layered signatures, bulk deconvolution algorithms have emerged as scalable alternatives to single-cell RNA-seq, though they remain sensitive to single-cell reference choices, data scaling (linear vs log), and stromal purity bias @addala2024computational. Here, we present a multi-layered transcriptomic and genomic analysis of 9 immunotherapy clinical trials from the iAtlas platform @iatlas2018 (n = 1,097) and 5 untreated indication-matched cohorts from The Cancer Genome Atlas (TCGA) @tcga2018 (n = 2,932). We integrate:
1. High-resolution single-cell deconvolution via `instaprism` @instaprism2024 to compare trial TME compositions with primary untreated baselines.
2. Machine learning classifiers predicting response from somatic mutations and identifying specific gene alterations (such as _TGM6_).
3. Quantitative subclonal dynamics modeling using Variant Allele Frequencies (VAF) to optimize TMB thresholds and predict TMB reliability.
4. Direct concordance benchmarking against continuous single-cell graph differential abundance (Milo) and multi-replicate perturbation stress-testing (@fig-study-overview).

#figure(
  image("../figures/dataset_overview/figure1_study_cohorts_overview.png", width: 100%),
  caption: [
    *Figure 1: Comprehensive multi-omic immunotherapy dataset compendium, single-cell references, and analytical deconvolution framework.*
    *(a)* Multi-cohort dataset ecosystem comprising a 40,002-cell pan-cancer scRNA-seq reference (6 major lineages, 22 cell types, 58 fine-grained sub-clusters) alongside the matched clinical validation cohort (Sade-Feldman et al., $n=51$ pre/post ICB melanoma biopsies, 16,291 CD45+ cells), 9 clinical trial cohorts ($n = 1,097$ patients) across melanoma, bladder, renal, pancreatic, and breast cancers, and 5 untreated primary tumor baseline cohorts from TCGA ($n = 2,932$ bulk tumors) defining $K=100$ unperturbed microenvironmental compositional centroids.
    *(b)* Multi-omic profiling layers spanning single-cell and bulk transcriptomics, whole-exome somatic mutations, subclonal variant allele frequency (VAF) spectra, tumor mutational burden (TMB), predictive response drivers (such as _TGM6_ alterations), standardized RECIST v1.1 clinical response criteria, and longitudinal treatment timing.
    *(c)* Unified analytical architecture executing single-cell continuous neighborhood testing (Milo kNN graph GLM), high-resolution pseudobulk deconvolution (`instaprism`), an automated multi-replicate distortion simulation suite ($N=5$, 485 runs across 8 single-variable modes, 4 compound regimes, and 2D interaction surfaces), and subclonal VAF joint optimization.
    *(d)* Key biomarker and methodological synthesis, illustrating empirical clinical concordance in matched melanoma ($rho = 0.853$), identification of concordant responder (naive B cells) and non-responder (macrophages) biomarkers alongside collinear cytotoxic T-cell discordance, primary deconvolution failure modes (patient marker dysregulation and cell size asymmetry), and quantification of the metastatic composition generalization gap.
  ]
) <fig-study-overview>

