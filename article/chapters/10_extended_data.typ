= Extended Data & Supplementary Information

#v(1em)

#figure(
  image("../figures/single-cell/qc_metadata_distributions.png", width: 90%),
  caption: [Supplementary Figure S1: Single-cell subsampling metadata and QC distributions. Distributions of total unique molecular identifier (UMI) counts, detected genes, mitochondrial percentage, and ribosomal protein percentage across the 40,002 subsampled cells.]
) <fig-supp-s1>

#pagebreak()

#figure(
  image("../figures/single-cell/umap_granularities.png", width: 95%),
  caption: [Supplementary Figure S2: UMAP projection of Leiden clustering resolutions. Evaluation across clustering resolutions (0.1 to 1.5) on the reference dataset showcasing lineage and subpopulation segregation.]
) <fig-supp-s2>

#pagebreak()

#figure(
  image("../figures/single-cell/umap_tumor_inference.png", width: 95%),
  caption: [Supplementary Figure S3: UMAP projection of dual-method tumor cell inference. Overlap of epithelial marker gene expression (left) and inferred chromosomal copy number variation (CNV) scores (right) used to identify malignant cells.]
) <fig-supp-s3>

#pagebreak()

#figure(
  image("../figures/single-cell/subclustering_metrics_distributions.png", width: 95%),
  caption: [Supplementary Figure S4: Daniela Witten Selective Inference p-values. Truncated normal selective inference p-values computed for all sub-clustering configurations, demonstrating statistically significant cluster separation ($p < 10^(-16)$).]
) <fig-supp-s4>

#pagebreak()

#figure(
  image("../figures/deconvolution/deconv_distributions_kmeans_subcluster_res_0.5.png", width: 100%),
  caption: [Supplementary Figure S5: Trial versus TCGA baseline deconvolution distributions. Side-by-side horizontal boxplots comparing the 58 sub-cluster fractions between the 9 trial cohorts and matched untreated TCGA indication cohorts.]
) <fig-supp-s5>

#pagebreak()

#figure(
  image("../figures/deconvolution/deconv_split_response_kmeans_subcluster_res_0.5.png", width: 100%),
  caption: [Supplementary Figure S6: Cell fraction distributions split by clinical response. Boxplot grid displaying K-means 0.5 sub-cluster fractions split by RECIST clinical response (Responders, green vs Non-Responders, red).]
) <fig-supp-s6>

#pagebreak()

#figure(
  image("../figures/deconvolution/deconv_split_tmb_kmeans_subcluster_res_0.5.png", width: 100%),
  caption: [Supplementary Figure S7: Cell fraction distributions split by tumor mutational burden. Boxplot grid displaying sub-cluster fractions split by cohort median DNA-TMB (High TMB, blue vs Low TMB, orange).]
) <fig-supp-s7>

#pagebreak()

#figure(
  image("../figures/deconvolution/kmeans_cluster_sizes.png", width: 95%),
  caption: [Supplementary Figure S8: TCGA K-means ($K=100$) cluster sizes and dominant cancer types. Sample count distributions across the 100 K-means centroids colored by primary cancer indication.]
) <fig-supp-s8>

#pagebreak()

#figure(
  image("../figures/deconvolution/kmeans_centroids_heatmap.png", width: 95%),
  caption: [Supplementary Figure S9: TCGA K-means centroid profiles. Heatmap of the 100 reference centroids across the 58 sub-clusters characterizing primary baseline architectures.]
) <fig-supp-s9>

#pagebreak()

#figure(
  image("../figures/deconvolution/iatlas_classification_umap.png", width: 95%),
  caption: [Supplementary Figure S10: UMAP embedding of nearest-centroid cancer type classification. Projection of the 1,097 trial samples onto TCGA reference centroids, illustrating the generalization gap in metastatic samples.]
) <fig-supp-s10>

#pagebreak()

#figure(
  image("../figures/deconvolution/deconv_prediction_metrics_grid.png", width: 95%),
  caption: [Supplementary Figure S11: Deconvolution response prediction metrics across feature selection strategies. Comparison of ROC AUC, PR AUC, and MCC across raw features, ANOVA univariate selection, and PCA dimensionality reduction.]
) <fig-supp-s11>

#pagebreak()

#figure(
  grid(
    columns: (1fr, 1fr),
    gutter: 1.2em,
    align(center)[
      #image("../figures/mutations/overfitting_underfitting_rosenberg.png", width: 100%)
      #v(0.3em)
      #text(size: 8pt)[(A) Rosenberg Bladder Cohort]
    ],
    align(center)[
      #image("../figures/mutations/overfitting_underfitting_melanoma.png", width: 100%)
      #v(0.3em)
      #text(size: 8pt)[(B) Combined Melanoma Cohort]
    ]
  ),
  caption: [Supplementary Figure S12: Model complexity and cross-validation performance. Comparison of training versus test set ROC AUC across Random Forest, Adaline, and MLP classifiers evaluating genomic overfitting.]
) <fig-supp-s12>

#pagebreak()

#figure(
  grid(
    columns: (1fr, 1fr),
    gutter: 1.2em,
    align(center)[
      #image("../figures/vaf-tmb/vaf_density_valley_detection.png", width: 100%)
      #v(0.3em)
      #text(size: 8pt)[(A) KDE Valley Detection]
    ],
    align(center)[
      #image("../figures/vaf-tmb/estimated_vs_true_thresholds.png", width: 100%)
      #v(0.3em)
      #text(size: 8pt)[(B) Estimated vs True Optimal VAF]
    ]
  ),
  caption: [Supplementary Figure S13: Unsupervised VAF threshold determination. (A) Kernel Density Estimation valley detection identifying clonal boundaries. (B) Scatter comparison of unsupervised estimated thresholds against empirically optimized cutoffs.]
) <fig-supp-s13>

#pagebreak()

#figure(
  grid(
    columns: (1fr, 1fr),
    gutter: 1.2em,
    align(center)[
      #image("../figures/vaf-tmb/reliability_score_vs_performance.png", width: 100%)
      #v(0.3em)
      #text(size: 8pt)[(A) TRS vs Predictive Strength]
    ],
    align(center)[
      #image("../figures/vaf-tmb/dna_vs_rna_tmb_performance.png", width: 100%)
      #v(0.3em)
      #text(size: 8pt)[(B) DNA-TMB vs Expressed RNA-TMB]
    ]
  ),
  caption: [Supplementary Figure S14: TMB reliability and RNA-TMB performance. (A) Positive correlation (Pearson $r = 0.527$) between TMB Reliability Score (TRS) and predictive strength. (B) Cross-validated ROC AUC comparison between DNA-TMB and Expressed RNA-TMB.]
) <fig-supp-s14>

#pagebreak()

#figure(
  grid(
    columns: (1fr, 1fr),
    gutter: 1.2em,
    align(center)[
      #image("../figures/vaf-tmb/odds_ratios_forest_plot_MCC.png", width: 100%)
      #v(0.3em)
      #text(size: 8pt)[(A) Forest Plot (MCC Objective)]
    ],
    align(center)[
      #image("../figures/vaf-tmb/odds_ratios_forest_plot_F1.png", width: 100%)
      #v(0.3em)
      #text(size: 8pt)[(B) Forest Plot (F1 Objective)]
    ]
  ),
  caption: [Supplementary Figure S15: Forest plots of Odds Ratios for optimized TMB/VAF splits. Fisher's Exact Test Odds Ratios and 95% Confidence Intervals for cutoffs optimized under (A) MCC and (B) F1-Score objectives.]
) <fig-supp-s15>

#pagebreak()

#figure(
  grid(
    columns: (1fr, 1fr),
    gutter: 1.2em,
    align(center)[
      #image("../figures/signatures/signature_distributions.png", width: 100%)
      #v(0.3em)
      #text(size: 8pt)[(A) Signature Score Distributions]
    ],
    align(center)[
      #image("../figures/signatures/signature_overlap_heatmap.png", width: 100%)
      #v(0.3em)
      #text(size: 8pt)[(B) Jaccard Overlap Heatmap]
    ]
  ),
  caption: [Supplementary Figure S16: ImmunoCompass transcriptional signatures. (A) Density distributions of signature scores across trial cohorts. (B) Jaccard similarity heatmap evaluating overlap and redundancy among ImmunoCompass gene sets.]
) <fig-supp-s16>

#pagebreak()

#figure(
  image("../figures/milo-deconv/perturbation_sensitivity_curves.png", width: 95%),
  caption: [Supplementary Figure S17: Multi-replicate single perturbation sensitivity curves ($N=5$ replicates per point). (A) Spearman rank correlation degradation (Mean ± SD). (B) Directional sign concordance (%) with dashed coin-flip baseline. (C) Ground truth deconvolution fidelity ($rho_("true","deconv")$). (D) Maximum degradation ranking sorted by peak correlation loss across 8 single distortion modes.]
) <fig-supp-s17>

#pagebreak()

#figure(
  image("../figures/milo-deconv/synthetic_perturbation_scatter_grid.png", width: 95%),
  caption: [Supplementary Figure S18: Sade-Feldman style scatter comparison grid for single perturbations. Eight panels displaying Baseline versus 6 distinct single distortion modes at representative moderate strength ($I approx 0.60$), colored by four-quadrant diagnostic classification with linear regression trends and correlation metrics.]
) <fig-supp-s18>
