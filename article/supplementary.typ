#set page(
  paper: "a4",
  margin: (x: 2cm, y: 2.5cm),
  header: align(right)[#text(size: 8pt, fill: gray)[Supplementary Materials | TME Deconvolution & Genomics]],
  footer: [
    #align(center)[#text(size: 8pt)[S#counter(page).display("1")]]
  ]
)

#set text(
  font: "Liberation Serif",
  size: 10pt,
  fill: rgb("#2C3E50")
)

#set par(justify: true, leading: 0.6em)

#align(center)[
  #text(size: 18pt, weight: "bold", fill: rgb("#1F3A60"))[
    Supplementary Materials
  ]
  #v(0.5em)
  #text(size: 12pt, weight: "medium")[
    Integrating Single-Cell Deconvolution, Somatic Mutations, and Subclonal Dynamics to Predict Immunotherapy Response across Clinical Trials and Baseline Cohorts
  ]
  #v(1em)
  #text(size: 10.5pt, style: "italic")[
    Antigravity AI Coding Assistant & Lu Jo Hae Team
  ]
  #v(0.5em)
  #text(size: 9pt, fill: gray)[
    Lu Jo Hae Lab, Department of Advanced Agentic Coding, Google DeepMind
  ]
]

#v(2.5em)

#show heading: it => [
  #set text(fill: rgb("#1F3A60"))
  #v(1.5em)
  #it
  #v(0.8em)
]

= Section S1: Single-Cell Reference Exploration & QC

#v(1em)
#align(center)[
  #image("../output/single-cell-exploration/qc_metadata_distributions.svg", width: 85%)
  #v(0.5em)
  #block(width: 85%)[
    #set text(size: 9pt)
    #set par(justify: true)
    *Supplementary Figure S1: Single-cell subsampling metadata and QC distributions.* Showing the distributions of total unique molecular identifier (UMI) counts, detected genes, mitochondrial gene percentage, and ribo-protein percentage across the subsampled cells.
  ]
]

#pagebreak()

= Section S2: Leiden Clustering Resolutions

#v(1em)
#align(center)[
  #image("../output/single-cell-exploration/umap_granularities.svg", width: 95%)
  #v(0.5em)
  #block(width: 90%)[
    #set text(size: 9pt)
    #set par(justify: true)
    *Supplementary Figure S2: UMAP projection of Leiden clustering resolutions.* Contrast between different Leiden resolution configurations (from 0.1 to 1.5) on the subsampled single-cell reference, showcasing the progression of cell partitioning.
  ]
]

#pagebreak()

= Section S3: Tumor Cell Detection

#v(1em)
#align(center)[
  #image("../output/single-cell-exploration/umap_tumor_inference.svg", width: 95%)
  #v(0.5em)
  #block(width: 90%)[
    #set text(size: 9pt)
    #set par(justify: true)
    *Supplementary Figure S3: UMAP projection of dual-method tumor cell inference.* Shows the overlap of (left) epithelial marker expression and (right) chromosomal copy number variation (CNV) proxy scores used to identify malignancy.
  ]
]

#pagebreak()

= Section S4: Selective Inference & Sub-Clustering

#v(1em)
#align(center)[
  #image("../output/single-cell-exploration/subclustering_metrics_distributions.svg", width: 95%)
  #v(0.5em)
  #block(width: 90%)[
    #set text(size: 9pt)
    #set par(justify: true)
    *Supplementary Figure S4: Daniela Witten Selective Inference p-values.* Truncated normal selective p-values computed for all sub-clustering configurations, showing highly significant separation ($p < 10^{-16}$).
  ]
]

#pagebreak()

= Section S5: Trial vs. Baseline K-Means Sub-Clusters

#v(1em)
#align(center)[
  #image("../output/single-cell-exploration/deconv_distributions_kmeans_subcluster_res_0.5.svg", width: 100%)
  #v(0.5em)
  #block(width: 95%)[
    #set text(size: 9pt)
    #set par(justify: true)
    *Supplementary Figure S5: Trial vs TCGA baseline deconvolution distributions for K-means 0.5 sub-clusters.* Comprehensive side-by-side horizontal boxplots comparing the 58 sub-cluster fractions between the 9 trial cohorts and their matched untreated TCGA indication comparators.
  ]
]

#pagebreak()

= Section S6: Cell Fractions Split by Response (K-Means 0.5)

#v(1em)
#align(center)[
  #image("../output/single-cell-exploration/deconv_split_response_kmeans_subcluster_res_0.5.svg", width: 100%)
  #v(0.5em)
  #block(width: 95%)[
    #set text(size: 9pt)
    #set par(justify: true)
    *Supplementary Figure S6: Cell fraction distributions split by immunotherapy response.* Boxplot grid displaying K-means 0.5 sub-cluster fractions split by clinical response (Responders, green vs Non-Responders, red).
  ]
]

#pagebreak()

= Section S7: Cell Fractions Split by TMB (K-Means 0.5)

#v(1em)
#align(center)[
  #image("../output/single-cell-exploration/deconv_split_tmb_kmeans_subcluster_res_0.5.svg", width: 100%)
  #v(0.5em)
  #block(width: 95%)[
    #set text(size: 9pt)
    #set par(justify: true)
    *Supplementary Figure S7: Cell fraction distributions split by tumor mutational burden (TMB).* Boxplot grid displaying K-means 0.5 sub-cluster fractions split by cohort median DNA-TMB (High TMB, blue vs Low TMB, orange).
  ]
]

#pagebreak()

= Section S8: Unsupervised Composition Clustering

#v(1em)
#align(center)[
  #image("../output/single-cell-exploration/kmeans_cluster_sizes.svg", width: 95%)
  #v(0.5em)
  #block(width: 90%)[
    #set text(size: 9pt)
    #set par(justify: true)
    *Supplementary Figure S8: TCGA K-means ($K=100$) cluster sizes and dominant cancer types.* Bar plot showing the distribution of sample counts across the 100 K-means centroids, with bars colored by dominant cancer type.
  ]
]

#pagebreak()

= Section S9: Centroids Heatmap

#v(1em)
#align(center)[
  #image("../output/single-cell-exploration/kmeans_centroids_heatmap.svg", height: 85%)
  #v(0.5em)
  #block(width: 90%)[
    #set text(size: 9pt)
    #set par(justify: true)
    *Supplementary Figure S9: Heatmap of K-means cluster centroids.* Centroid fractions (y-axis: 100 clusters, sorted by dominant cancer type; x-axis: 22 cell type fractions) showing distinct compositional profiles.
  ]
]

#pagebreak()

= Section S10: iAtlas Cancer Type Classification

#v(1em)
#align(center)[
  #image("../output/single-cell-exploration/iatlas_classification_umap.svg", width: 95%)
  #v(0.5em)
  #block(width: 90%)[
    #set text(size: 9pt)
    #set par(justify: true)
    *Supplementary Figure S10: iAtlas actual vs predicted cancer type UMAP.* UMAP projection of the iAtlas trial samples colored by (left) actual trial cancer type and (right) K-means predicted cancer type, showing the generalization gap.
  ]
]

#pagebreak()

= Section S11: Baseline Response Prediction Metrics

#v(1em)
#align(center)[
  #image("../output/single-cell-exploration/deconv_prediction_metrics_grid.svg", width: 95%)
  #v(0.5em)
  #block(width: 90%)[
    #set text(size: 9pt)
    #set par(justify: true)
    *Supplementary Figure S11: 6-panel grid of baseline response prediction metrics.* Grouped barplots comparing PR AUC, Accuracy, Precision, Recall, F1, and MCC against the random baseline (gray) for all trial cohorts.
  ]
]

#pagebreak()

= Section S12: Somatic Mutation Overfitting Plots

#v(1em)
#grid(
  columns: (1fr, 1fr),
  gutter: 1em,
  align(center)[
    #image("../output/iAtlas-response-mutations/Rosenberg-iAtlas/overfitting_underfitting_comparison.svg", width: 95%)
    #v(0.2em)
    #text(size: 8pt)[(A) Rosenberg Bladder]
  ],
  align(center)[
    #image("../output/iAtlas-response-mutations/Combined-Melanoma/overfitting_underfitting_comparison.svg", width: 95%)
    #v(0.2em)
    #text(size: 8pt)[(B) Combined Melanoma]
  ]
)
#v(0.5em)
#block(width: 100%)[
  #set text(size: 9pt)
  #set par(justify: true)
  *Supplementary Figure S12: Somatic mutation response prediction training vs validation curves.* Overlay of training AUC (red) and validation AUC (blue) across 23 feature selection configurations for Random Forest, Adaline, and MLP-4 classifiers. Shows severe overfitting for neural networks on large feature sizes (e.g. all genes) but stable fits on restricted features (e.g. fish_10).
]

#pagebreak()

= Section S13: Subclonal VAF Valleys & Threshold Comparisons

#v(1em)
#grid(
  columns: (1.1fr, 0.9fr),
  gutter: 1em,
  align(center)[
    #image("../output/iAtlas-vaf-tmb-reliability/vaf_density_valley_detection.svg", width: 95%)
    #v(0.2em)
    #text(size: 8pt)[(A) VAF Density Valleys]
  ],
  align(center)[
    #image("../output/iAtlas-vaf-tmb-reliability/estimated_vs_true_thresholds.svg", width: 95%)
    #v(0.2em)
    #text(size: 8pt)[(B) Unsupervised vs. True Optimal VAF]
  ]
)
#v(0.5em)
#block(width: 100%)[
  #set text(size: 9pt)
  #set par(justify: true)
  *Supplementary Figure S13: Subclonal VAF valleys and threshold comparison.* (A) Unsupervised KDE valley detection identifying subclonal/clonal transition boundaries across cohorts. (B) Scatter comparison showing estimated unsupervised VAF thresholds (KDE, GMM, Otsu) plotted against true optimal VAF thresholds maximizing clinical response association.
]

#pagebreak()

= Section S14: TMB Reliability Score (TRS) & RNA-TMB

#v(1em)
#grid(
  columns: (1fr, 1fr),
  gutter: 1em,
  align(center)[
    #image("../output/iAtlas-vaf-tmb-reliability/reliability_score_vs_performance.svg", width: 95%)
    #v(0.2em)
    #text(size: 8pt)[(A) TRS vs TMB Predictive Strength]
  ],
  align(center)[
    #image("../output/iAtlas-vaf-tmb-reliability/dna_vs_rna_tmb_performance.svg", width: 95%)
    #v(0.2em)
    #text(size: 8pt)[(B) DNA-TMB vs. Expressed RNA-TMB]
  ]
)
#v(0.5em)
#block(width: 100%)[
  #set text(size: 9pt)
  #set par(justify: true)
  *Supplementary Figure S14: TMB reliability and RNA-TMB performance.* (A) Positive association (Pearson $r = 0.527$) between TMB Reliability Score (TRS) and TMB predictive strength ($|"AUC" - 0.5|$). (B) Out-of-fold cross-validated ROC AUC comparison showing that filtering somatic mutations by bulk expression (Expressed RNA-TMB, green) degrades response prediction compared to standard DNA-TMB (blue).
]

#pagebreak()

= Section S15: Odds Ratios Forest Plots (Optimized TMB/VAF splits)

#v(1em)
#grid(
  columns: (1fr, 1fr),
  gutter: 1em,
  align(center)[
    #image("../output/iAtlas-vaf-tmb-optimization/odds_ratios_forest_plot_MCC.svg", width: 95%)
    #v(0.2em)
    #text(size: 8pt)[(A) Forest Plot (MCC Optimization)]
  ],
  align(center)[
    #image("../output/iAtlas-vaf-tmb-optimization/odds_ratios_forest_plot_F1.svg", width: 95%)
    #v(0.2em)
    #text(size: 8pt)[(B) Forest Plot (F1 Optimization)]
  ]
)
#v(0.5em)
#block(width: 100%)[
  #set text(size: 9pt)
  #set par(justify: true)
  *Supplementary Figure S15: Forest plots of Odds Ratios for optimized TMB/VAF splits.* Showing Fisher's Exact Test Odds Ratios and their 95% Confidence Intervals for splits optimized under (A) Matthews Correlation Coefficient (MCC) and (B) F1-Score objectives. Significant clinical associations ($p < 0.05, "OR" > 1.0$) are highlighted in green.
]

#pagebreak()

= Section S16: ImmunoCompass Signatures

#v(1em)
#grid(
  columns: (1fr, 1fr),
  gutter: 1em,
  align(center)[
    #image("../output/immuno-compass-signatures/signature_distributions.svg", width: 95%)
    #v(0.2em)
    #text(size: 8pt)[(A) Signature Score Distributions]
  ],
  align(center)[
    #image("../output/immuno-compass-signatures/signature_overlap_heatmap.svg", width: 95%)
    #v(0.2em)
    #text(size: 8pt)[(B) Jaccard Overlap Heatmap]
  ]
)
#v(0.5em)
#block(width: 100%)[
  #set text(size: 9pt)
  #set par(justify: true)
  *Supplementary Figure S16: ImmunoCompass transcriptional signatures.* (A) Density distributions of signature scores across clinical cohorts. (B) Jaccard index heatmap evaluating overlap and redundancy among the ImmunoCompass signature gene sets.
]
