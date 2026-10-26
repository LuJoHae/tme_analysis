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
To assess when pseudobulk deconvolution faithfully recapitulates single-cell differential abundance, we formalized the mathematical mapping between continuous graph testing and constrained linear simplex inversion @dann2022differential, @chu2022cell (@fig-how-milopy-works).

#figure(
  placement: auto,
  image("../figures/deconvolution/how_milopy_works.png", width: 100%),
  caption: [Methodological architecture of Milo / milopy differential abundance testing on single-cell k-nearest neighbor graphs. (1) kNN graph construction in PCA latent space and index vertex sampling with medoid refinement, bypassing discrete clustering boundaries. (2) Overlapping neighborhood definition, sparse binary incidence matrix $bold(M) in {0, 1}^(C times V)$, pairwise graph distance, and cross-cohort neighborhood composition $bold(C) = bold(B)^T bold(A)$. (3) Aggregation of single cells into patient-level count table $bold(N) in bb(N)_0^(V times J)$ with TMM library size normalization, eliminating single-cell pseudo-replication through hierarchical patient replicate modeling. (4) Quasi-likelihood Negative Binomial GLM with Empirical Bayes dispersion shrinkage, Spatial FDR p-value weighting via local graph connectivity density, and dual outputs (volcano plot and single-cell UMAP $log_2"FC"$ projection).]
) <fig-how-milopy-works>

In single-cell continuous neighborhood testing (Milo), an index graph $G = (V, E)$ is constructed in latent principal component space ($k = 30$ nearest neighbors). For each sampled index vertex $v$, cell counts across patient samples $j in {1, ..., J}$ are modeled using a Negative Binomial generalized linear model with log link:
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

== 2.8 Hierarchical Deconvolution & Reference Collinearity Resolution (BayesPrism Framework)
When fine phenotypic subsets share marker programs (reference collinearity, $bold(phi)_(s_2) = bold(phi)_(s_1) + bold(delta)$ with $norm(bold(delta))_2 << norm(bold(phi)_(s_1))_2$), the reference Gram matrix $bold(Phi)^T bold(Phi)$ becomes ill-conditioned ($kappa(bold(Phi)) -> infinity$). Evaluating the Fisher information matrix $cal(I)(bold(theta))$ along the normalized sibling difference vector $bold(v) = 1/sqrt(2) (bold(e)_(s_1) - bold(e)_(s_2))$ reveals that $bold(v)^T cal(I)(bold(theta)) bold(v) = O(norm(bold(delta))_2^2) -> 0$. By the Cramér-Rao bound, the posterior covariance explodes along $bold(v)$:
$ "Var"(hat(theta)_(s_1) - hat(theta)_(s_2)) >= 2 bold(v)^T cal(I)^(-1) bold(v) = O(1 / norm(bold(delta))_2^2) -> infinity $
inducing strong negative cross-talk: $"Cov"(hat(theta)_(s_1), hat(theta)_(s_2)) approx -sqrt("Var"(hat(theta)_(s_1)) "Var"(hat(theta)_(s_2)))$.

To circumvent this singularity while preserving intra-lineage transcriptomic heterogeneity, BayesPrism @chu2022cell defines a surjective partition mapping $pi: cal(S) -> cal(T)$ from fine cell states $cal(S) = {1, ..., S}$ into broad cell types $cal(T) = {1, ..., T}$ ($T << S$), formalized by the linear aggregation operator $bold(M) in {0, 1}^(T times S)$ where $M_(t,s) = 1$ if $pi(s) = t$ and $0$ otherwise. 

Under this framework, deconvolution executes across four coordinated stages:
1. *State-Level Latent Sampling*: Unobserved count matrices $bold(Z)_(n,g,s)$ and state fractions $bold(theta)_(n,s)$ are sampled from a fine-state Dirichlet-Multinomial Gibbs sampler ($bold(Phi)^("state") in bb(R)^(S times G)$), capturing patient-specific state mixtures without forcing an artificial single centroid.
2. *Null-Space Invariance & Variance Cancellation*: State estimates are marginalized to broad types:
$ Z_(n,g,t) = sum_(s in cal(S)_t) Z_(n,g,s) = (bold(M) bold(Z)_(n,g,dot))_t, quad theta_(n,t)^((0)) = sum_(s in cal(S)_t) theta_(n,s) = (bold(M) bold(theta)_n)_t $
Crucially, the collinear difference direction lies identically in the null space of the aggregation operator:
$ bold(M) bold(v) = 1/sqrt(2) (bold(M) bold(e)_(s_1) - bold(M) bold(e)_(s_2)) = 1/sqrt(2) (bold(e)_t - bold(e)_t) = bold(0) ==> bold(v) in "Null"(bold(M)) $
Consequently, the divergent variance terms of order $O(1/norm(bold(delta))_2^2)$ cancel against the negative covariance terms:
$ "Var"(theta_(n,t)^((0))) = sum_(s in cal(S)_t) "Var"(hat(theta)_(n,s)) + 2 sum_(s < s' in cal(S)_t) "Cov"(hat(theta)_(n,s), hat(theta)_(n,s')) = O(1) $
3. *Broad Reference Refinement*: Marginalized latent counts are pooled across the cohort ($Z_(g,t) = sum_n Z_(n,g,t)$) to update broad type profiles $bold(psi)_t = "Softmax"(log bold(phi)_t^("type") + bold(gamma)_t)$ under a zero-mean Gaussian shrinkage prior $gamma_(t,g) ~ cal(N)(0, sigma^2)$. MAP optimization prevents overfitting and platform bias without inflating collinear state degrees of freedom.
4. *Final Non-Collinear Deconvolution*: Final fractions $bold(theta)_n^((f))$ are sampled strictly on the well-conditioned simplex $Delta^(T-1)$ ($kappa(bold(Psi)) << 100$), eliminating all collinear instability.

== 2.9 Collinearity-Aware Regularized Deconvolution on the Probability Simplex
While BayesPrism's broad-type aggregation operator $bold(M)$ successfully cancels the null-space singularity at the lineage level ($bold(v) in "Null"(bold(M))$), it fundamentally collapses fine-grained cell states into coarse clusters, relinquishing intra-lineage phenotypic resolution. To retain state-level discrimination without suffering condition number collapse, we formulated a collinearity-aware regularized objective defined natively on the canonical probability simplex $Delta^(S-1) = {bold(theta) in bb(R)_+^S : sum_s theta_s = 1}$:
$ min_(bold(theta) in Delta^(S-1)) cal(L)(bold(theta)) = D_("KL")(bold(x), bold(Phi)^T bold(theta)) + lambda_("lap")/2 bold(theta)^T bold(L) bold(theta) + lambda_("fuse") sum_(i < j) W_(i j) psi_epsilon (theta_i - theta_j) $
where $D_("KL")(bold(x), bold(mu)) = sum_g (mu_g - x_g log mu_g)$ represents the unnormalized Kullback-Leibler deviance under Poisson counting statistics, with expected expression $bold(mu) = bold(Phi)^T bold(theta)$.

The regularizer incorporates two transcriptomically informed geometric operators:
1. *Graph-Laplacian Manifold Energy*: Let $bold(W) in bb(R)_+^(S times S)$ denote the symmetric cell-state transcriptomic affinity matrix, with entries $W_(i j) = ("Corr"(bold(phi)_i, bold(phi)_j)_+)^2$ for distinct states ($i eq.not j$) sharing lineage assignment ($M_(t, i) = M_(t, j) = 1$) and $0$ otherwise. The unnormalized graph Laplacian $bold(L) = bold(D) - bold(W)$ (where $D_(i i) = sum_j W_(i j)$) induces a Dirichlet energy:
$ bold(theta)^T bold(L) bold(theta) = 1/2 sum_(i, j) W_(i j) (theta_i - theta_j)^2 $
which penalizes divergent assignments along collinear directions where $W_(i j) -> 1$.
2. *Smoothed Fused Lasso*: To prevent over-smoothing across genuinely distinct states, the pseudo-Huber smoothed $ell_1$ penalty $psi_epsilon(u) = sqrt(u^2 + epsilon^2)$ (with smoothing parameter $epsilon = 10^(-5)$) encourages grouping of closely related states without non-differentiable singularity at $u=0$.

To solve this constrained non-linear program efficiently without manifold distortion, we developed a projected Fast Iterative Shrinkage-Thresholding Algorithm (FISTA) equipped with Armijo backtracking line search @beck2009fista. In each iteration $k$, the candidate state is evaluated via:
$ bold(theta)^((k+1)) = Pi_(Delta^(S-1)) (bold(y)^((k)) - alpha_k nabla cal(L)(bold(y)^((k)))) $
where $Pi_(Delta^(S-1))$ denotes the exact Euclidean projection onto the probability simplex computed in $O(S log S)$ time via the sorting algorithm of Wang and Carreira-Perpiñán @wang2013projection.

To guarantee that regularization strictly prevents condition number explosion without manual hyperparameter tuning or destructive overregularization, we derived an automated *Adaptive Deficit Spectral Calibration*. In well-conditioned or moderately collinear regimes (e.g., $kappa(bold(H)_"nom") approx 2 times 10^3$), first-order simplex solvers easily resolve likelihood curvature; imposing an excessively aggressive uniform target (such as $kappa = 10.0$) forces the Graph Laplacian penalty to dwarf the data deviance, inducing midpoint shrinkage bias where unequal true proportions are artificially pulled toward equality ($theta_1 approx theta_2$). To eliminate this overregularization floor, our adaptive calibration injects curvature strictly proportional to the eigenvalue deficit beyond an operational target ceiling $kappa_("target")$ (default $2,000.0$):
$ sigma_("target") = (sigma_max (bold(H)_"nom")) / kappa_("target"), quad Delta sigma = max(0, sigma_("target") - sigma_min (bold(H)_"nom")), quad lambda_("lap") = (Delta sigma) / overline(d) $
where $bold(H)_"nom" = bold(Phi) "diag"(1 / bold(mu)_0) bold(Phi)^T$ is the nominal Poisson Hessian evaluated at the barycentric prior $bold(theta)_0 = 1/S bold(1)$, and $overline(d) = 1/S "Tr"(bold(L))$ is the average node degree. Because the Hessian along the collinear null vector $bold(v) = 1/sqrt(2)(bold(e)_i - bold(e)_j)$ satisfies:
$ bold(v)^T nabla^2 cal(L)(bold(theta)) bold(v) = bold(v)^T bold(H)_"nom" bold(v) + lambda_("lap") bold(v)^T bold(L) bold(v) = O(norm(bold(delta))_2^2) + 2 lambda_("lap") W_(i j) >= 2 lambda_("lap") W_(i j) $
the regularizer injects strictly positive curvature $+2 lambda_("lap") W_(i j)$ orthogonal to the likelihood manifold, guaranteeing that the effective condition number satisfies $kappa(nabla^2 cal(L)) <= kappa_("target")$ when states become identical ($norm(bold(delta))_2 -> 0$), while vanishing ($lambda_("lap") approx 0$) when nominal conditioning is already sufficient.


== 2.10 Finite-Sample Controlled Stability Selection & Regularization Budgeting
To reliably extract transcriptomic markers predictive of immunotherapy response without succumbing to high-dimensional sample variance and false discovery inflation ($p >> n$), we implemented stability selection @meinshausen2010stability equipped with complementary pairs subsampling @shah2013variable.

Given feature matrix $bold(X) in bb(R)^(n times p)$ and binary response vector $bold(y) in {0, 1}^n$, the data are repeatedly partitioned into complementary subsamples of size $floor(n/2)$ across $B$ bootstrap cycles ($2B$ subsamples indexed by $b in {1, ..., 2B}$). For any model selection algorithm parameterized by regularization vector $bold(lambda) in Lambda$, let $hat(S)^bold(lambda)(I_b) subset.eq {1, ..., p}$ denote the set of active features selected on subsample $I_b$. The empirical selection probability for feature $k$ is given by:
$ hat(Pi)_k^bold(lambda) = 1/(2B) sum_(b=1)^(2B) bb(I)(k in hat(S)^bold(lambda)(I_b)) $
The stable feature set is defined by thresholding the maximum selection probability across the regularization path $Lambda$:
$ hat(S)^("stable") = {k in {1, ..., p} : max_(bold(lambda) in Lambda) hat(Pi)_k^bold(lambda) >= pi_("thr")} $
where $pi_("thr") in (0.5, 1.0]$ represents the stability selection threshold (set to $0.60$ for single cohorts and $0.70$ for multi-cohort pooled analyses).

Under the exchangeability assumption for noise variables, the Per-Family Error Rate (PFER, defined as the expected number of falsely selected noise variables $bb(E)[V]$) is bounded by:
$ bb(E)[V] <= 1/(2 pi_("thr") - 1) q_Lambda^2 / p $
where $q_Lambda = bb(E)[|hat(S)^bold(lambda)(I)|]$ represents the average number of features selected by the base procedure across the path, and $p$ is the total feature space ($p = 822$ across 5 immune-oncology gene sets). Under unimodality or $r$-concavity conditions on the selection probability distributions, Shah and Samworth @shah2013variable established tighter error bounds:
$ bb(E)[V] <= (q_Lambda^2 / p) C(pi_("thr"), theta) $
where $C(pi_("thr"), theta) < 1 / (2 pi_("thr") - 1)$ for $pi_("thr") > 0.5$.

To strictly prevent over-selection and algorithmic runaway along decreasing regularization paths, we implemented a rigorous *Regularization Budgeting Engine*. The theoretical active set capacity $q_("target")$ is calibrated from user-specified nominal PFER bounds:
$ q_("budget") = floor(sqrt("PFER"_("target") dot p dot (2 pi_("thr") - 1))) $
During path evaluation across ordered regularization penalties $lambda_1 > lambda_2 > ... > lambda_L$, the empirical budget expenditure $overline(q)(lambda_l) = 1/(2B) sum_b |hat(S)^(lambda_l)(I_b)|$ is tracked dynamically. Iteration terminates immediately once $overline(q)(lambda_l) >= q_("budget")$, ensuring that the finite-sample error guarantees are preserved unconditionally across all fitters.

== 2.11 Multi-Cohort & Multi-Paradigm Fitter Framework
Heterogeneous clinical cohorts introduce systematic technical batch effects, institutional sequencing protocols, and distinct patient histology baselines. To evaluate marker stability across diverse model inductive biases and multi-cohort confounding structures, we engineered a unified, modular fitter framework encompassing 13 distinct feature selection algorithms across five distinct mathematical paradigms:

1. *Convex Regularized GLMs*:
  - *Lasso* ($ell_1$ penalized logistic regression @tibshirani1996regression): Coordinate descent minimizing binomial deviance with $ell_1$ penalty, enforcing feature-level sparsity.
  - *Elastic Net* ($ell_1 + ell_2$ mixture @zou2005regularization): Balances sparse selection with group shrinkage through penalty $alpha norm(bold(beta))_1 + ((1 - alpha)/2) norm(bold(beta))_2^2$ ($alpha = 0.5$), retaining correlated gene clusters.

2. *Non-Convex & Adaptive Penalization*:
  - *SCAD* (Smoothly Clipped Absolute Deviation) & *MCP* (Minimax Concave Penalty): Local linear approximations (LLA) mitigating asymptotic estimation bias of $ell_1$ regularization for large coefficients.
  - *Adaptive Lasso*: Two-stage reweighted $ell_1$ penalization with weight vector $w_j = 1/|hat(beta)_j^("init")|^gamma$ ($gamma = 1.0$), satisfying oracle selection properties.
  - *Hard Thresholding*: Iterative parameter projection retaining strictly top-$k$ gradient magnitude entries.

3. *Gradient Tree Ensembles*:
  - *LightGBM & XGBoost*: Histogram-based gradient boosted decision trees utilizing split gain and feature importance coverage to evaluate non-linear synergistic interactions.
  - *Random Forest*: Bagged ensemble of decorrelated classification trees using Gini impurity reduction.

4. *Cohort-Aware Confounder Partialling & Hierarchical Mixed Models*:
  - *Frisch-Waugh-Lovell Partialled Lasso (`LassoFWL`)*: To purge inter-cohort batch effects without inflating gene degrees of freedom, we apply the Frisch-Waugh-Lovell theorem @frisch1933partial @lovell1963seasonal. Let $bold(C) in bb(R)^(n times K)$ denote the one-hot cohort indicator matrix. Defining the orthogonal cohort projection operator $bold(M)_bold(C) = bold(I)_n - bold(C)(bold(C)^T bold(C))^(-1)bold(C)^T$, the residualized expression matrix $bold(X)^* = bold(M)_bold(C) bold(X)$ and response vector $bold(y)^* = bold(M)_bold(C) bold(y)$ reside strictly in the orthogonal complement of the cohort subspace, purging study-level mean offsets prior to $ell_1$ selection.
  - *Mixed Effects Random Forest (`MERF`)*: Decomposes patient outcome into a shared non-linear fixed-effects transcriptome function $f(bold(x))$ and cohort-specific random intercepts $bold(b) ~ cal(N)(0, sigma_b^2 bold(I)_K)$ @hajjem2014mixed:
  $ y_(i j) = f(bold(x)_(i j)) + b_j + epsilon_(i j), quad i in {1, ..., n_j}, quad j in {1, ..., K} $
  Parameters are estimated via Expectation-Maximization (EM), alternating between forest residual fitting and empirical Bayes random effect updating.

5. *Ordered Weighted $ell_1$ Regularization*:
  - *OSCAR* & *SLOPE*, detailed below.

== 2.12 Ordered Weighted $ell_1$ Norms: OSCAR Octagonal Clustering & SLOPE FDR Control
High-dimensional biomarker discovery in transcriptomics is severely impaired by multi-gene co-expression modules and false discovery inflation. To jointly resolve gene clustering and exact false discovery rate control, we implemented Ordered Weighted $ell_1$ (OWL) regularized objectives:
$ min_(bold(beta) in bb(R)^p) cal(L)(bold(beta)) = 1/(2 n) norm(bold(y) - bold(X) bold(beta))_2^2 + Omega_bold(w)(bold(beta)) $
where the ordered weighted norm $Omega_bold(w)(bold(beta))$ is defined for non-negative non-increasing weights $w_1 >= w_2 >= ... >= w_p >= 0$ by:
$ Omega_bold(w)(bold(beta)) = sum_(j=1)^p w_j |beta|_((j)) $
where $|beta|_((1)) >= |beta|_((2)) >= ... >= |beta|_((p))$ represents the sequence of absolute coefficients sorted in descending order.

1. *OSCAR (Octagonal Shrinkage and Clustering Algorithm for Regression)*:
Bondell and Reich @bondell2008simultaneous showed that simultaneous sparsity and pairwise clustering ($|beta_i - beta_j| -> 0$) can be represented as an OWL norm with affine linear weights:
$ w_j^("OSCAR") = lambda_1 + lambda_2 (p - j), quad j in {1, ..., p} $
When coefficients are unequal, the differential penalty creates octagonal polytope facets that force identical magnitude and sign on correlated covariates, forming data-driven gene clusters.

2. *SLOPE (Sorted L-One Penalized Estimation)*:
Bogdan et al. @bogdan2015slope established that setting the OWL weights according to the decaying Gaussian quantile sequence:
$ w_j^("SLOPE") = lambda dot Phi^(-1)(1 - (q_("FDR") dot j)/(2 p)), quad j in {1, ..., p} $
(where $Phi^(-1)$ is the standard normal quantile and $q_("FDR") in (0, 1)$ is the target false discovery rate, set to $0.10$) guarantees finite-sample False Discovery Rate (FDR) control at level $q_("FDR")$ for independent or orthogonal designs, while adapting to unknown coefficient sparsity.

*Fast Proximal Gradient Solver via PAVA*:
Evaluating the proximal operator for $Omega_bold(w)(bold(v)) = "argmin"_(bold(u)) { 1/2 norm(bold(u) - bold(v))_2^2 + Omega_bold(w)(bold(u)) }$ requires projecting the sorted absolute values onto the monotone non-increasing cone. Let $bold(P) in {0, 1}^(p times p)$ be the permutation matrix sorting $|bold(v)|$ such that $bold(z) = bold(P)|bold(v)|$ satisfies $z_1 >= z_2 >= ... >= z_p >= 0$. The proximal solution satisfies:
$ "prox"_(Omega_bold(w))(bold(v)) = "diag"("sgn"(bold(v))) bold(P)^T "PAVA"_+(bold(z) - bold(w)) $
where $"PAVA"_+$ applies the Pool Adjacent Violators Algorithm (isotonic regression) to enforce monotonicity $theta_1 >= theta_2 >= ... >= theta_p >= 0$, resolving pairwise clusters in $O(p log p)$ time. Optimization is accelerated via the Fast Iterative Shrinkage-Thresholding Algorithm (FISTA) with backtracking line search @beck2009fista.


== 2.13 Unified `tme_datasets` Registry and Multi-Stage Computational Provenance
To eliminate cohort-specific preprocessing drift, manual path coupling, and gene-identifier mismatch across single-cell and bulk transcriptomic modalities, all clinical cohorts, reference atlases, and immuno-oncology gene panels were centralized in the declarative `packages/tme_datasets` registry (@fig-tme-datasets-pipeline-dag). Every dataset is governed by an immutable `DatasetSpec` record specifying modality (`SINGLE_CELL` or `BULK_RNA`), canonical accession and study aliases, primary organ and disease indication, sequencing chemistry, quality-control bounds (`QualityControlSpec`), and harmonized clinical endpoints. Across all registered cohorts, gene identifiers are mapped to canonical HGNC protein-coding symbols under Ensembl release 111, and clinical outcomes are standardized onto binary RECIST v1.1 response indicators ($Y_n = 1$ for Complete or Partial Response, `CR`/`PR`, or Major Pathological Response in neoadjuvant cohorts; $Y_n = 0$ for Stable or Progressive Disease, `SD`/`PD`).

#figure(
  placement: auto,
  image("../figures/dataset_overview/figure_tme_datasets_pipeline_dag.png", width: 100%),
  caption: [
    *Unified `tme_datasets` registry queries and multi-stage computational provenance directed acyclic graph (DAG).*
    *(a) Tier I — Declarative Registry, Providers & Gene Panels (`packages/tme_datasets`):* Three foundational query domains expose (i) 9 single-cell ICB response cohorts (`SC-9`, $3.60 times 10^6$ cells across melanoma, pan-cancer, colorectal, triple-negative breast, non-small cell lung, head-and-neck squamous, and pancreatic ductal cancers) via `load_dataset()`, `resolve_milopy_cohorts()`, and `load_backed()`; (ii) 11 bulk RNA-seq clinical trials plus 4 indication/pan-cancer pooled strata (`BULK-15`, $N = 1,058$ evaluable biopsies across 9 iAtlas trials and 2 NSCLC cohorts) via `load_dataset()` and `load_iatlas_cohort_or_combined()`; and (iii) 17+ auxiliary single-cell reference atlases, 5 untreated TCGA primary tumor baselines ($N = 2,932$), and curated gene collections (`IMMUNOTHERAPY_GENE_PANEL`, `IMMUNE_CHECKPOINT_GENES`, Bagaev TME signatures) via `load_geneset_collection()` and `list_registered_datasets()`.
    *(b) Tier II — Primary Query, QC & `.h5ad` Serialization Scripts (`scripts/manuscript_pipeline/`):* `1a_download_and_preprocess_all.py`, `1b_batch_preprocess_cohorts.py`, and `1c_run_qc_pipeline_all_single_cell.py` `[1A–1C]` execute Ensembl v111 symbol reconciliation and write atomic sparse `AnnData` caches (`data/preprocessed/{dataset_id}.h5ad`) alongside per-script summaries under `output/manuscript_pipeline/<script_stem>/`; in parallel, `2a1_sample_and_analyze_cohorts.py` and `2a2_combine_cohorts_and_analyze.py` `[2A]` query `SC-9` to compute single-cell kNN Negative Binomial GLM differential abundance (`output/manuscript_pipeline/2a1_sample_and_analyze_cohorts/*.parquet`), while `2b_run_stability_selection_iatlas.py` `[2B]` queries `BULK-15` and the 101-gene IO panel across 12 cohorts $times$ 13 fitters (`output/manuscript_pipeline/2b_run_stability_selection_iatlas/*.parquet`).
    *(c) Tier III — Chained Analytical Compute Engines:* Downstream scripts in `scripts/manuscript_pipeline/` consume upstream `.h5ad` and `.parquet` artifacts deterministically and isolate outputs in `output/manuscript_pipeline/<script_stem>/`: `3a_sade_feldman_pseudobulk_deconv_benchmark.py` `[3A]` pairs the 9 single-cell `.h5ad` caches `[1A]` with `cell_level_scores.parquet` `[2A]` to construct 355 Leiden ($r = 2.0$) reference centroids; `3b_pan_cohort_unified_deconv_benchmark.py` `[3B]` Ward-clusters those 355 states into the 42-state `Consensus-Meta42` pan-cancer atlas (`M01`–`M42`); `3c_iatlas_instaprism_milopy_benchmark.py` `[3C]` integrates the 11 bulk `.h5ad` caches `[1A]` with references from `[3A]` and `[3B]` to execute bulk linearization, InstaPrism EM deconvolution across 13 references $times$ 15 cohorts, 5-family ML response prediction, zero-shot `LOCO` transfer, DerSimonian–Laird `Pan-iAtlas-Meta` concordance, 168-gene marker specificity scoring, and $8 times 8$ / $42 times 42$ non-linear epistasis; while `3d_curate_expert_gene_list.py` `[3D]`, `3e1_benchmark_iatlas_icb_predictors.py`, and `3e2_run_hpo_pipeline.py` `[3E]` synthesize stability dossiers and HPO benchmarks.
    *(d) Tier IV — Publication Figure, Table & Manuscript Synthesis:* `4b_plot_iatlas_instaprism_milopy_figures.py` `[4B]` consumes all 11 Parquet tables emitted by `[3C]` to render the 11 cross-modality benchmark and full-atlas figures (`output/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures/` and `article/figures/milo-deconv/iatlas_*.{svg,png}`), joined by `[4A]` (`4a1`–`4a7`), `[4C]` (`4c1`–`4c5`), and `[4D]` (`4d_generate_tme_datasets_pipeline_dag.py`) in `article/Makefile` for 300 DPI `rsvg-convert` rasterization and Typst compilation into `article/main.pdf`.
  ]
) <fig-tme-datasets-pipeline-dag>

As mapped in @fig-tme-datasets-pipeline-dag, all 26 manuscript pipeline scripts are topologically ordered inside `scripts/manuscript_pipeline/` (`1a` through `4d`), and every script writes its artifacts into a dedicated `output/manuscript_pipeline/<script_stem>/` subdirectory across four strictly decoupled tiers:
1. *Tier I Registry & Providers*: Exposes typed query primitives over 9 single-cell ICB cohorts (`GSE120575` @sadefeldman2018defining, `CELLxGENE_7b20c613` @gondal2025curated, `CELLxGENE_05a8c945` @marteau2026singlecell, `CELLxGENE_6f9de485` @yan2026singlecell, `GSE207422` @hu2023singlecell, `GSE243013` @liu2025multiomic, `GSE233203` @hong2025singlecell, `GSE200996` @luoma2022molecular, and `GSE316195` @raufi2026singlenucleus; totaling $3.60 times 10^6$ cells), 9 bulk RNA-seq trials from CRI iAtlas (`Hugo` @hugo2016genomic, `Riaz` @riaz2017tumor, `Liu` @liu2019association, `Gide` @gide2019clinical, `Rosenberg` @rosenberg2016atezolizumab, `Padron` @padron2022pancreatic, `Anders` @anders2022evaluating, `McDermott` @mcdermott2018clinical, and `Choueiri` @choueiri2016pembrolizumab; $N = 1,015$ evaluable pre-treatment patients), 2 independent NSCLC anti-PD-(L)1 bulk cohorts (`Jung` / `GSE135222` @jung2019dna, $n = 27$; and `Cho` / `GSE126044` @cho2020genome, $n = 16$), 4 multi-trial indication/pan-cancer pools (`Melanoma-Combined`, $n = 338$; `RenalCell-Combined`, $n = 263$; `NSCLC-Combined`, $n = 43$; `PanCancer-Combined`, $n = 1,015$), 5 TCGA baselines ($N = 2,932$), and curated IO gene panels.
2. *Tier II Ingestion, Single-Cell Milo DA & Stability Selection*: `1a_download_and_preprocess_all.py` `[1A]` materializes standardized sparse `AnnData` caches in `data/preprocessed/`. `2a1_sample_and_analyze_cohorts.py` `[2A]` queries `resolve_milopy_cohorts()` to fit patient-aware kNN Negative Binomial GLMs, exporting neighborhood and cell-level response log-fold-changes (`output/manuscript_pipeline/2a1_sample_and_analyze_cohorts/*.parquet`), while `2b_run_stability_selection_iatlas.py` `[2B]` exports finite-sample budgeted stability trajectories (`output/manuscript_pipeline/2b_run_stability_selection_iatlas/*.parquet`).
3. *Tier III Chained Reference Construction & Cross-Modality Benchmark*: `3a_sade_feldman_pseudobulk_deconv_benchmark.py` `[3A]` ingests the 9 `.h5ad` caches `[1A]` and `milopy` cell scores `[2A]` to derive 355 fine-grained Leiden ($r = 2.0$) centroids; `3b_pan_cohort_unified_deconv_benchmark.py` `[3B]` aggregates these centroids into the 42-state `Consensus-Meta42` reference (`output/manuscript_pipeline/3b_pan_cohort_unified_deconv_benchmark/*.parquet`); and `3c_iatlas_instaprism_milopy_benchmark.py` `[3C]` consumes `[1A]`, `[3A]`, and `[3B]` to produce 11 cross-modality Parquet tables in `output/manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark/`. Parallel modules `3d_curate_expert_gene_list.py` `[3D]` and `3e1_benchmark_iatlas_icb_predictors.py` / `3e2_run_hpo_pipeline.py` `[3E]` consume `[2B]` and `tme_datasets` bulk cohorts to emit biomarker dossiers and predictor benchmarks.
4. *Tier IV Declarative Visualization & Manuscript Build*: Pure plotting generators (`4b_plot_iatlas_instaprism_milopy_figures.py` `[4B]`, `4a2_generate_combined_figures.py` / `4a7_assemble_immunotherapy_2x2_figure.py` `[4A]`, `4c2_combine_all_stability_paths.py` `[4C]`, and `4d_generate_tme_datasets_pipeline_dag.py` `[4D]`) read strictly from Tier II–III Parquet outputs without re-querying raw data, exporting layered SVGs into `output/manuscript_pipeline/<script_stem>/` rendered at 300 DPI via `rsvg-convert` in `article/Makefile`.


== 2.14 Hierarchical Single-Cell Reference Tiers, Bulk Linearization, and Cross-Modality ICB Prediction
To evaluate how single-cell reference granularity, indication matching, and collinearity consolidation govern bulk RNA-seq deconvolution and clinical response prediction across independent trials, we constructed a four-tier hierarchy of 13 single-cell reference matrices $bold(Phi) in bb(R)_+^(G times K)$ from the 9 single-cell ICB cohorts ($3.60 times 10^6$ cells):
1. *Tier 1A — Cohort-Specific Fine-State References (9 matrices, $K in [29, 54]$)*: Within each single-cell cohort $c in {1, ..., 9}$, Harmony-integrated PCA neighborhoods were partitioned via Leiden clustering at resolution $r = 2.0$, yielding 355 total fine-grained states (`GSE120575`: $K = 38$; `CELLxGENE_7b20c613`: $K = 48$; `CELLxGENE_05a8c945`: $K = 39$; `CELLxGENE_6f9de485`: $K = 42$; `GSE207422`: $K = 38$; `GSE243013`: $K = 33$; `GSE233203`: $K = 34$; `GSE200996`: $K = 29$; `GSE316195`: $K = 54$). For each state $k$, we recorded both its mean linear count-per-ten-thousand profile $bold(phi)_k$, mean total UMI cell-size factor $S_k$, and patient-level `milopy` differential abundance ground truth $overline("logFC")_k$.
2. *Tier 1B — Indication-Unified References (2 matrices)*: Concatenation of shared-gene centroids across melanoma/pan-cancer single-cell cohorts (`Melanoma-PanCancer-Unified`, $K = 86$ states combining `GSE120575` and `CELLxGENE_7b20c613`) and across the three lung cancer cohorts (`NSCLC-Unified`, $K = 105$ states combining `GSE207422`, `GSE243013`, and `GSE233203`).
3. *Tier 2 — Uncollapsed Pan-Cancer Concatenation (`Unified-All9`, $K = 355$)*: Direct horizontal union $bold(Phi)^("All9") = [bold(Phi)^((1)) | dots.c | bold(Phi)^((9))]$ across all 355 cohort-specific states over the intersection of expressed protein-coding genes, intentionally retaining cross-study duplicate lineages to quantify unregularized multi-atlas collinearity.
4. *Tier 3 — Ward-Clustered Pan-Cancer Consensus Atlas (`Consensus-Meta42`, $K = 42$)*: Pairwise cosine distances $D_(a b) = 1 - " cos"(tilde(bold(phi))_a, tilde(bold(phi))_b)$ were computed across all 355 $z$-standardized log-centroids over highly variable genes and partitioned into $M = 42$ consensus meta-states (`M01`–`M42`) via Ward's minimum-variance linkage. Each meta-state signature $bold(phi)_m^("meta")$ and single-cell ground-truth effect $overline("logFC")_m^("meta")$ were aggregated via softmax cosine-affinity weights ($tau = 12.0$) weighted by constituent cell counts, collapsing cross-cohort collinear duplicates onto a well-conditioned pan-cancer basis while preserving specialized niches (`CXCL13+PDCD1+ Tex`, `TCF7+IL7R+ Stem CD8`, `CXCL9+STAT1+ M1-TAM`, `MS4A1+CD79A+ TLS B`, `JCHAIN+MZB1+ Plasma`, `ISLR+COL5A2+ myCAF`, and `FAP+TGFB1+ Mat. CAF`).

*Bulk Scale Detection, Exact Linearization, and Whole-Tumor Frame Alignment*:
Because linear and Multinomial deconvolution likelihoods operate strictly in additive RNA abundance space ($bold(b)_n approx bold(Phi) bold(theta)_n^("mRNA")$), feeding log-transformed bulk expression directly into deconvolution induces severe Jensen's inequality distortion. Prior to deconvolution, each bulk cohort matrix $bold(X) in bb(R)^(N times G)$ was audited automatically: cohorts stored on a logarithmic scale ($max_(n, g) X_(n, g) < 50$, including the iAtlas log2-normalized expression matrices and log-FPKM GEO series) were inverted to linear scale via $L_(n, g) = max(2^(X_(n, g)) - 1, 0)$, whereas raw linear count matrices were retained in linear space, followed by sample-wise Transcripts-Per-Million simplex normalization:
$ "TPM"_(n, g) = (10^6 L_(n, g)) / (sum_(g'=1)^G L_(n, g')) $
Furthermore, because `GSE120575` @sadefeldman2018defining profiles FACS-sorted $"CD45"^+$ leukocytes whereas unfractionated clinical bulk RNA-seq samples whole tumor biopsies containing malignant parenchyma, endothelial vessels, and cancer-associated fibroblasts, immune-only references were augmented with pan-cancer non-immune stromal/malignant background components during bulk EM estimation (`InstaPrism` @instaprism2024, 100 EM iterations) and subsequently renormalized onto the evaluated state simplex after dividing mRNA fractions by cell-size factors $S_k$:
$ hat(theta)_(n, k) = (hat(f)_(n, k)^("mRNA") / S_k) / (sum_(m=1)^K hat(f)_(n, m)^("mRNA") / S_m) $

*Compositional Feature Engineering, Multi-Model Supervised Prediction, and Zero-Shot LOCO Transfer*:
To map simplex-constrained proportions $hat(bold(theta))_n in Delta^(K-1)$ into Euclidean feature space without closure bias, inferred state fractions were transformed via Aitchison's centered log-ratio ($"CLR"_(n, k) = log((hat(theta)_(n, k) + delta) / g(hat(bold(theta))_n))$ with geometric mean $g(dot)$ and pseudocount $delta = 10^(-4)$) alongside variance-stabilizing arcsine-square-root coordinates $arcsin(sqrt(hat(theta)_(n, k)))$, concatenated and $z$-standardized strictly inside training folds. Across all $13 times 15 = 195$ (reference $times$ bulk cohort/pool) pairs, we benchmarked five supervised classifier families @scikit-learn:
1. $ell_2$-penalized Logistic Regression (`ridge_l2`),
2. Elastic-Net Logistic Regression (`elastic_net`, $ell_1$ ratio $= 0.5$),
3. $ell_1$-penalized Lasso Logistic Regression (`lasso_l1`),
4. Ridge-regularized Linear Probability Model (`linear_reg`), and
5. Non-linear Random Forest (`random_forest`, 300 trees, balanced subsample class weights, minimum leaf size 3).

Within-cohort discrimination was quantified via $5 times 5$ repeated stratified cross-validation (25 out-of-fold evaluations per cell; adaptive 3-fold CV when minority class count $< 5$), reporting Receiver Operating Characteristic Area Under the Curve (`ROC-AUC`) with Hanley–McNeil 95% confidence intervals, Precision–Recall AUC (`PR-AUC`), and prevalence-calibrated excess precision gain $Delta "PR-AUC" = "PR-AUC" - pi_("pos")$, where $pi_("pos") = N_("resp") / N$ is the trial responder prevalence. To test prospective cross-trial generalizability without cohort-specific refitting, we executed zero-shot *Leave-One-Cohort-Out (`LOCO`)* transfer across all 11 individual clinical trials: each model was trained on the remaining 10 trials (encoding cancer indication as one-hot covariates during training) and evaluated zero-shot on the held-out trial's compositional profiles.

*Random-Effects Meta-Analytic Cross-Modality Concordance and Marker Specificity*:
For each bulk trial $b$ and reference state $k$, univariate responder association $hat(beta)_("deconv", b, k)$ and its asymptotic variance $v_(b, k)$ were estimated via standardized point-biserial logit transformation (Section 2.5). Across the 9 independent iAtlas trials ($N = 1,015$), pan-cancer random-effects meta-analytic effect sizes (`Pan-iAtlas-Meta`) were synthesized via the DerSimonian–Laird estimator @dersimonian1986meta:
$ hat(beta)_("meta", k) = (sum_(b=1)^9 w_(b, k)^* hat(beta)_("deconv", b, k)) / (sum_(b=1)^9 w_(b, k)^*), quad w_(b, k)^* = (v_(b, k) + hat(tau)_k^2)^(-1) $
where $hat(tau)_k^2 = max(0, (Q_k - (B - 1)) / C_k)$ captures inter-trial heterogeneity. Cross-modality concordance between bulk $hat(beta)_("deconv", k)$ and single-cell `milopy` $overline("logFC")_k$ was assessed via Pearson linear correlation $r$, Spearman rank correlation $rho$, directional sign concordance among active states, and exact two-sided label-permutation tests ($B = 2,000$ permutations). For biological interpretability across both the focused 14-state poles and the complete 42-state `Consensus-Meta42` atlas (`168` state-gene pairs), each state's top-4 protein-coding HGNC marker genes (excluding ribosomal `RPS`/`RPL`, mitochondrial `MT-`, and immunoglobulin variable `IGKV`/`IGHV`/`IGLV` artifacts) were ranked by one-vs-rest reference specificity:
$ log_2 "FC"_("spec", k, g) = log_2 ((phi_(k, g) + epsilon_0) / (1/(K - 1) sum_(m eq.not k) phi_(m, g) + epsilon_0)) $
and cross-checked against unconfounded bulk expression responder Welch $t$-statistics across `PanCancer-Combined` ($N = 1,015$).


== 2.15 Non-Linear Pairwise Epistasis, Random Forest 2D Partial Dependence, and Conditional Stromal Gating
While linear classifiers capture additive shifts in cell-state log-ratios, anti-tumor immunity requires spatially coordinated cellular partnerships that can be physically or biochemically gated by immunosuppressive stroma. To systematically identify cell states that act primarily as conditional context modifiers rather than univariate main-effect markers, we regressed pan-cancer out-of-fold Random Forest Gini importance $I_k^("RF")$ against absolute single-cell abundance shift $|overline("logFC")_k|$ across all $K = 42$ `Consensus-Meta42` states, flagging positive-residual outlier states with high non-linear split leverage despite modest marginal effect sizes.

To dissect pairwise synergy and antagonism both among the top outlier/anchor states ($8 times 8$ focused submatrix) and exhaustively across the entire `Consensus-Meta42` ecosystem ($binom(42, 2) = 861$ unique pairs; $42 times 42 = 1,764$ ordered entries), we combined two complementary statistical metrics across the pooled `PanCancer-Combined` cohort ($N = 1,015$):
1. *Trial-Adjusted Multiplicative Logistic Epistasis ($z_(i j)$)*: For every unordered pair of $z$-standardized CLR state abundances $(z_(n, i), z_(n, j))$, we fitted a cohort-stratified logistic interaction model controlling for trial baseline log-odds via one-hot trial fixed effects $bold(c)_n in {0, 1}^8$:
$ "logit" Pr(Y_n = 1 | z_(n, i), z_(n, j), bold(c)_n) = beta_0 + beta_i z_(n, i) + beta_j z_(n, j) + gamma_(i j) (z_(n, i) z_(n, j)) + bold(alpha)^T bold(c)_n $
recording the two-sided Wald interaction statistic $z_(i j) = hat(gamma)_(i j) / "SE"(hat(gamma)_(i j))$. Positive interactions ($z_(i j) > +2$) indicate supra-additive co-activation synergy, whereas negative interactions ($z_(i j) < -2$) indicate conditional stromal or myeloid gating (where accumulation of state $j$ attenuates the responder slope of state $i$).
2. *Friedman–Popescu 2D Partial-Dependence Interaction Strength ($H_(i j)$)*: To quantify model-free non-additive probability curvature learned by tree ensembles @friedman2008predictive, let $hat(F)_(i j)(z_i, z_j)$ denote the centered two-variable partial dependence function of the pan-cancer Random Forest and $hat(F)_i(z_i), hat(F)_j(z_j)$ denote its centered marginal partial dependence functions evaluated over quantile grids. Friedman's interaction statistic measures the fraction of joint partial-dependence variance inexplicable by additive main effects:
$ H_(i j)^2 = (sum_(n=1)^N [hat(F)_(i j)(z_(n, i), z_(n, j)) - hat(F)_i(z_(n, i)) - hat(F)_j(z_(n, j))]^2) / (sum_(n=1)^N hat(F)_(i j)^2(z_(n, i), z_(n, j))) $
Global epistatic hub centrality for each state $i in {1, ..., 42}$ was summarized by its mean absolute interaction magnitude across all 41 partners, $overline(|z|)_i = 1/41 sum_(j eq.not i) |z_(i j)|$, paired with mean $overline(H)_i = 1/41 sum_(j eq.not i) H_(i j)$.

Finally, to translate continuous epistasis into clinically interpretable patient stratification regimes across all $N = 1,015$ patients, we constructed $3 times 3$ conditional tertile gating tables (`T1: Low`, `T2: Moderate`, `T3: High`) at two biological scales: (i) a *focused lymphoid-humoral effector score* (`M04 + M30 + M29`) stratified by a *focused fibrotic/M2 stromal barrier score* (`M17 myCAF + M32 M2-TAM`), and (ii) a *pan-atlas composite Responder TME score* (mean $z$-standardized abundance across all 26 responder-leaning states with $hat(beta)_("meta", k) > 0$) stratified by a *pan-atlas composite Non-Responder TME barrier score* (mean $z$-standardized abundance across all 16 non-responder-leaning states with $hat(beta)_("meta", k) < 0$). Within each cell of the $3 times 3$ grid, empirical responder rates $hat(p)_(r, c)$ were accompanied by exact Wilson score 95% confidence intervals and stratum-conditional point-biserial correlations $r(Y, S_("resp") | S_("barrier") in T_r)$.
