== 3.2 Trial vs. Baseline Deconvolution
Deconvoluting the bulk cohorts using `instaprism` revealed a strong, highly significant rank-order conservation of cell type abundances between trial cohorts and matching TCGA baseline indications (@tab-deconv-corr, @fig-supp-s5). Spearman rank correlation coefficients ($rho$) spanned between *0.72 and 0.86* ($p < 10^(-4)$) across metastatic melanoma (`Hugo`, `Riaz`), bladder urothelial carcinoma (`Rosenberg`), breast cancer (`Anders`), and clear cell renal cell carcinoma (`McDermott`).

#figure(
  table(
    columns: (1.2fr, 1fr, 1.2fr, 1.5fr),
    inset: 4.5pt,
    align: center,
    [*Trial Cohort*], [*TCGA Baseline*], [*Pearson $r$*], [*Spearman $rho$*],
    [Hugo], [SKCM], [0.571], [0.807 ($p < 10^(-5)$)],
    [Riaz], [SKCM], [0.509], [0.787 ($p < 10^(-4)$)],
    [Gide], [SKCM], [-0.057], [0.162 ($p = 0.47$)],
    [Rosenberg], [BLCA], [0.365], [0.723 ($p < 10^(-3)$)],
    [Anders], [BRCA], [0.461], [0.860 ($p < 10^(-6)$)],
    [McDermott], [KIRC], [0.470], [0.831 ($p < 10^(-5)$)]
  ),
  caption: [Instaprism deconvolution correlation between clinical trial cohorts and matching TCGA primary baselines.]
) <tab-deconv-corr>

Despite this high rank conservation, significant compositional divergence was observed. In particular, trial biopsies demonstrated a marked elevation of myeloid and cancer-associated fibroblast fractions. Strikingly, the `Gide-iAtlas` cohort (heavily pre-treated metastatic melanoma biopsies) exhibited near-zero correlation with `TCGA-SKCM` ($rho = 0.162, p = 0.47$), underscoring massive therapy-induced microenvironmental remodeling. Examining cell fractions split by clinical response (@fig-supp-s6) and DNA-TMB status (@fig-supp-s7) revealed distinct lineage-specific shifts.

== 3.3 Unsupervised Generalization Gap
Classifying the 1,097 trial samples against $K=100$ TCGA-derived reference centroids revealed a severe generalization gap (@tab-gen-gap), yielding an overall classification accuracy of only *26.34%*. 

#figure(
  table(
    columns: (1.4fr, 1fr, 1fr, 1fr),
    inset: 4.5pt,
    align: center,
    [*Cancer Type*], [*Precision*], [*Recall*], [*F1-Score*],
    [SKCM (Melanoma)], [0.89], [0.52], [0.66],
    [KIRC (Renal)], [0.46], [0.28], [0.34],
    [BRCA (Breast)], [0.04], [0.97], [0.08],
    [BLCA (Bladder)], [0.00], [0.00], [0.00],
    [PAAD (Pancreas)], [0.00], [0.00], [0.00]
  ),
  caption: [Performance of nearest-centroid cancer type classification of trial samples mapped onto TCGA baseline centroids.]
) <tab-gen-gap>

Bladder (`BLCA`) and pancreatic (`PAAD`) samples failed entirely to map to their corresponding primary centroids (recall = 0%), driven by metastatic site organotropism, stromal overgrowth, and pre-treatment chemotherapy signatures (@fig-supp-s8, @fig-supp-s9, @fig-supp-s10).
