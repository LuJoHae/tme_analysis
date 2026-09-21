= 2. Materials and Methods

== 2.1 Single-Cell Reference & Quality Control
We subsampled *40,002 cells* from a pan-cancer scRNA-seq reference across diverse tumor and stromal lineages. Quality control filtering evaluated unique molecular identifier (UMI) counts, detected genes per cell, and mitochondrial and ribosomal gene proportions (@fig-supp-s1). Leiden clustering was evaluated across resolution parameters ranging from 0.1 to 1.5 (@fig-supp-s2), with resolution 0.5 chosen to delineate 22 transcriptionally distinct cell types. Malignancy was assessed by integrating epithelial marker gene expression with inferred chromosomal copy number variations (CNV) (@fig-supp-s3).

Within each major lineage, fine-grained cell states (58 sub-clusters) were defined using Silhouette-optimized K-means clustering. Centroid separations were rigorously validated using truncated normal selective inference p-values according to the Daniela Witten framework (@fig-supp-s4) @selectiveinference2023. Parallel deconvolution of bulk RNA-seq data from 9 clinical trial cohorts and 5 TCGA baseline cohorts was executed using `instaprism` with 50 iterations @instaprism2024.

== 2.2 Unsupervised Composition Classifier & Generalization Gap
To assess microenvironmental similarity and divergence, we clustered 2,932 TCGA baseline deconvolution profiles into $K=100$ K-means clusters (@fig-supp-s8, @fig-supp-s9). Each cluster centroid was assigned a dominant cancer type based on majority vote. We then projected the 1,097 iAtlas clinical trial samples into this latent compositional space and classified each trial sample by finding its nearest TCGA centroid via Euclidean distance (@fig-supp-s10).

== 2.3 Somatic Mutations & Response Prediction
Binary somatic mutation indicators were extracted across all sequenced genes. We trained Random Forest, Adaline, and Multi-Layer Perceptron (MLP) classifiers @scikit-learn using Stratified 5-Fold Cross-Validation across 23 feature configurations to evaluate predictive power and quantify overfitting (@fig-supp-s11, @fig-supp-s12). For candidate response-associated genes (such as _TGM6_), two-sided Fisher's Exact Tests were conducted to calculate Odds Ratios and significance, while Mann-Whitney U tests evaluated differences in TMB distributions.

== 2.4 Subclonal VAF & TMB Optimization
Tumor Mutational Burden (TMB) was calculated across 10 Variant Allele Frequency (VAF) thresholds ranging from 0.01 to 0.50. Joint optimization of VAF ($t_("vaf")$) and TMB ($t_("tmb")$) cutoff thresholds was conducted under a regularized objective with a Log-Normal prior:
$ f_("reg") = f - lambda ((log(t_("vaf")) - log(0.05))^2 + (log(t_("tmb")) - log(10.0))^2) $
penalizing excessive deviations from established clinical baselines ($t_("vaf") = 0.05, t_("tmb") = 10.0$ mut/Mb) with penalty factor $lambda = 0.20$ (@fig-supp-s15).

Unsupervised VAF threshold estimators were evaluated using Kernel Density Estimation (KDE) valley detection, Otsu's binarization, and Gaussian Mixture Models (GMM) boundaries (@fig-supp-s13). To predict whether a given cohort's TMB profile yields reliable prognostic signal, we formulated the TMB Reliability Score (TRS):
$ "TRS" = "MER" times log_(10)("Median Depth") $
where MER denotes the somatic mutation expression rate in bulk RNA-seq (@fig-supp-s14). In parallel, published transcriptional signatures from the ImmunoCompass compendium were scored across cohorts (@fig-supp-s16).

== 2.5 Mathematical Modeling: Milo DA vs. Pseudobulk Deconvolution
To assess when pseudobulk deconvolution faithfully recapitulates single-cell differential abundance, we formalized the mathematical mapping between continuous graph testing and constrained linear simplex inversion @dann2022differential, @chu2022cell. 

In single-cell continuous neighborhood testing (Milo), an index graph $G = (V, E)$ is constructed in latent principal component space ($k = 15$ nearest neighbors). For each sampled index vertex $v$, cell counts across patient samples $j in {1, ..., J}$ are modeled using a Negative Binomial generalized linear model with log link:
$ log(bb(E)[N_(v,j)]) = beta_(0,v) + beta_v^("milo") Y_j + log("TotalCells"_j) $
where $Y_j in {0, 1}$ denotes binary clinical response, and $log("TotalCells"_j)$ functions as a library-size offset. For discrete cell cluster $k$, the cluster-level Milo effect size is summarized as $overline("logFC")_k = "median"_{v in C_k} (beta_v^("milo"))$.

Conversely, pseudobulk deconvolution models the bulk transcriptome $bold(b)_j in bb(R)^G$ as a convex combination of cell-type expression signatures $bold(Phi) in bb(R)^(G times K)$:
$ bold(b)_j approx bold(Phi) bold(f)_j^("mRNA"), quad bold(f)_j^("mRNA") in Delta^(K-1) $
Crucially, deconvolution inherently estimates mRNA mass fractions $f_(k,j)^("mRNA")$, which decouple from physical cell count fractions $p_(k,j)$ under variable per-cell mRNA content $S_k$:
$ f_(k,j)^("mRNA") = (p_(k,j) S_k) / (sum_(m=1)^K p_(m,j) S_m) $
To establish mathematical comparability with $overline("logFC")_k$, patient-level deconvolution fractions $hat(f)_(k,j)$ are standardized, and effect sizes $hat(beta)_k^("deconv")$ are estimated via standardized point-biserial logit transforms:
$ r_(p b, k) = "Corr"(Y_j, hat(f)_(k,j)), quad hat(beta)_k^("deconv") = (2 r_(p b, k)) / sqrt(1 - r_(p b, k)^2 + epsilon) $

== 2.6 Multi-Replicate Perturbation Suite & Compound Clinical Regimes
We developed an automated, modular simulation suite evaluating 8 distinct biological, technical, and transcriptional distortion modes across 5 independent replicates per condition ($N=5$, 485 total simulations):
1. *Cell Size / mRNA Asymmetry* ($S_k in [1times, 50times]$): Exponential mRNA yield disparity across lineages.
2. *Reference Collinearity* ($alpha in [0, 99%]$): Marker gene sharing between adjacent clusters, driving condition number $kappa(bold(Phi)) > 10^3$.
3. *Inter-Patient Shift* ($sigma in [0, 3.0]$): Log-normal patient-specific baseline expression drift.
4. *State Activation Confounding* ($gamma in [1times, 20times]$): Marker upregulation in responders without cell division.
5. *Unmodeled Ghost / Tumor Contamination* ($f_("ghost") in [0, 85%]$): Bulk reads originating from omitted parenchymal tumor.
6. *Sampling Sparsity* ($sigma in [0, 4.0]$): Biopsy-level cell yield inequality with severe Poisson dropouts.
7. *Ambient RNA ("Soup") Contamination* ($eta in [0, 40%]$): Droplet background RNA diffusion into single-cell profiles.
8. *Patient Marker Dysregulation* ($sigma_("dysreg") in [0, 2.0]$): Patient-specific repression of lineage markers in bulk tissue.

To model real-world clinical pathology where distortions co-occur, we parameterized 4 compound clinical regimes across intensity $alpha in [0.0, 1.0]$ ($120$ simulations):
- *Clinical Core Needle Biopsy*: Cell Size ($1.0 times alpha$) + Sparsity ($0.8 times alpha$) + Ghost Contamination ($0.7 times alpha$).
- *Inflamed Tumor Microenvironment*: Activation ($1.0 times alpha$) + Collinearity ($0.8 times alpha$) + Patient Shift ($0.6 times alpha$).
- *Single-Cell Technical Noise*: Ambient Soup ($1.0 times alpha$) + Sparsity ($0.6 times alpha$) + Patient Shift ($0.5 times alpha$).
- *Triple-Jeopardy Breakdown*: Cell Size ($1.0 times alpha$) + Activation ($1.0 times alpha$) + Ghost Contamination ($0.8 times alpha$).

In addition, a 2D Factorial Interaction sweep ($5 times 5$ grid of Cell Size $times$ Activation Confounding, $N=5$ reps, $125$ simulations) was executed to quantify non-linear interaction surfaces.

== 2.7 Matched Clinical Validation Pipeline (Sade-Feldman Cohort)
To empirically test our mathematical models against clinical reality, we benchmarked single-cell Milo DA directly against self-reference pseudobulk deconvolution in the Sade-Feldman metastatic melanoma cohort (GSE120575, $n=51$ biopsies, 12 immune cell states) @sadefeldman2018defining. Patient pseudobulks were reconstructed by aggregating raw counts per biopsy. Deconvolution was executed using `instaprism` with the single-cell cluster centroids as reference. For each cell state, directional concordance was classified into four diagnostic quadrants: Concordant Responder ($hat(beta)_k > 0.1, overline("logFC")_k > 0.1$), Concordant Non-Responder ($hat(beta)_k < -0.1, overline("logFC")_k < -0.1$), Discordant, or Concordant Neutral.
