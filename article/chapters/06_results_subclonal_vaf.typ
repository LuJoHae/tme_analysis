== 3.6 Subclonal VAF & TMB Threshold Optimization
Somatic mutation Variant Allele Frequencies (VAF) were significantly correlated with clinical response, yet the direction of association proved highly cohort-dependent: responders exhibited elevated VAF in Rosenberg bladder ($p = 3.3 times 10^(-291)$) but depressed VAF in McDermott renal cell carcinoma ($p = 1.9 times 10^(-12)$).

Joint optimization of the VAF calculation cutoff and TMB stratification threshold demonstrated sharp, cohort-specific association peaks (@fig-vaf-peaks, @tab-vaf-splits):
- `Rosenberg-iAtlas` (Bladder): Maximal response discrimination peaked at *VAF = 0.05* (Mann-Whitney U $p = 1.77 times 10^(-8)$, Cohen's $d = 0.693$, ROC AUC = 0.758).
- `Combined-Melanoma`: The association peaked at a clonal threshold of *VAF = 0.03* ($p = 0.0254$, Cohen's $d = 0.301$, ROC AUC = 0.591).
- `Padron-iAtlas` (Pancreatic): Optimal separation emerged at *VAF = 0.42* ($p = 0.0075$, Cohen's $d = -0.522$), where elevated clonal TMB paradoxically associated with non-response, consistent with immunosuppressive desmoplastic exclusion.

#figure(
  image("../figures/vaf-tmb/Rosenberg-iAtlas_vaf_tmb_association.png", width: 90%),
  caption: [Subclonal VAF and TMB response association landscape in the Rosenberg bladder cancer cohort, illustrating optimal separation at VAF = 0.05 and TMB = 10.0 mut/Mb.]
) <fig-vaf-peaks>

#figure(
  table(
    columns: (1.4fr, 1fr, 1fr, 1.6fr),
    inset: 4.5pt,
    align: center,
    [*Clinical Cohort*], [*Opt VAF*], [*Opt TMB*], [*Fisher OR (95% CI)*],
    [Rosenberg (BLCA)], [0.05], [10.0], [5.08 (3.01 - 8.79)],
    [Hugo (SKCM)], [0.26], [22.0], [10.43 (1.70 - 75.33)],
    [Riaz (SKCM)], [0.30], [21.0], [3.97 (1.28 - 14.15)]
  ),
  caption: [Regularized optimal VAF and TMB cutoffs maximizing odds of clinical response across immunotherapy cohorts.]
) <tab-vaf-splits>

Unsupervised *Kernel Density Estimation (KDE) valley detection* successfully identified optimal VAF thresholds without clinical outcome supervision (@fig-supp-s13), matching true optimal values closely in Rosenberg (estimated 0.02 vs empirical 0.05) and Anders (estimated 0.04 vs empirical 0.12).

Furthermore, our *TMB Reliability Score (TRS)* successfully anticipated cohort-level TMB predictive fidelity (Pearson *$r = 0.527$*, @fig-supp-s14 #text[(A)]). Evaluating *Expressed RNA-TMB* demonstrated that filtering somatic mutations by bulk expression failed to improve predictive value (ROC AUC remained flat or decreased, @fig-supp-s14 #text[(B)]), confirming that bulk RNA-seq lacks sufficient sensitivity to capture low-frequency immunogenic neoantigens presented on tumor MHC complexes. Systematic forest plots of optimized odds ratios across cohorts are provided in @fig-supp-s15.
