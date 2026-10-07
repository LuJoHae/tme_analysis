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
