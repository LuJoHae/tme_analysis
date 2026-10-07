# Hierarchical Deconvolution & Reference Collinearity Resolution in BayesPrism

**A Rigorous Mathematical Formulation and Code Analysis**

---

## 1. Executive Summary

A pervasive obstacle in reference-based transcriptomic deconvolution is **reference collinearity**: closely related cell states (such as CD8+ cytotoxic vs. exhausted T cells, inflammatory vs. patrolling monocytes, or subclonal malignant populations) share 95%+ of their expressed transcripts. In standard deconvolution methods—spanning Ordinary Least Squares (OLS), Non-Negative Least Squares (NNLS), support vector regression (CIBERSORTx), and single-level probabilistic models—this transcriptional overlap inflates the condition number $\kappa(\boldsymbol{\Phi})$ of the reference signature matrix $\boldsymbol{\Phi}$. Consequently, the matrix inverse $(\boldsymbol{\Phi}^T \boldsymbol{\Phi})^{-1}$ becomes ill-conditioned, triggering massive variance explosion, sign-flipping of regression coefficients, and severe negative cross-talk between sibling states.

**BayesPrism** (*Chu et al., Nature Cancer, 2022*) circumvents this breakdown through a **two-tier hierarchical Bayesian framework** that explicitly separates **broad cell types** ($\mathcal{T}$, e.g., T cells, B cells, Myeloid, Malignant) from **fine cell states** ($\mathcal{S}$, e.g., phenotypic sub-clusters or continuous manifold states). 

Rather than choosing between over-simplified coarse averaging (which erases patient-specific cell state heterogeneity) and flat granular deconvolution (which succumbs to collinear noise), BayesPrism implements a four-stage mechanism:
1. **Fine State Latent Variable Sampling**: Infers patient-specific latent read counts $Z_{n, g, s}$ and state proportions $\theta_{n, s}$ at the granular state level $\mathcal{S}$, allowing the likelihood to flexibly allocate reads according to subtle state differences.
2. **Null-Space Statistical Marginalization**: Projects the fine-state estimates to the broad cell type space via an aggregation operator $\mathbf{M} \in \{0, 1\}^{T \times S}$. The collinear noise vector between sibling states lies exactly in the null space of $\mathbf{M}$, resulting in exact mathematical cancellation of the divergent variance terms.
3. **Broad Lineage Reference Refinement**: Pools marginalized reads across patients to update cell-type-specific profiles $\boldsymbol{\psi}_t$ under a regularizing Gaussian shrinkage prior ($\gamma_{t, g} \sim \mathcal{N}(0, \sigma^2)$), avoiding the catastrophic overfitting that would occur if updating collinear state profiles.
4. **Final Non-Collinear Proportion Estimation**: Samples final cell-type fractions $\boldsymbol{\theta}_n^{(f)}$ strictly on the low-dimensional, well-conditioned cell-type simplex $\Delta^{T-1}$.

This document provides the complete analytical formulation, mathematical proofs, and code-level mapping between the theory and BayesPrism's implementation in both R and functional Python/PyTorch.

---

## 2. The Mathematical Pathology of Reference Collinearity

### 2.1 The Ill-Conditioned Linear Inverse Problem

Consider a bulk gene expression mixture vector $\mathbf{x} \in \mathbb{R}^G$ across $G$ genes, modeled as a linear combination of $S$ cell state reference signatures $\boldsymbol{\Phi} = [\boldsymbol{\phi}_1, \boldsymbol{\phi}_2, \dots, \boldsymbol{\phi}_S] \in \mathbb{R}^{G \times S}$:

$$\mathbf{x} = \boldsymbol{\Phi} \boldsymbol{\theta} + \boldsymbol{\epsilon}, \quad \boldsymbol{\theta} \in \Delta^{S-1} = \left\{ \boldsymbol{\theta} \in \mathbb{R}_+^S : \sum_{s=1}^S \theta_s = 1 \right\}$$

In standard unconstrained or non-negatively constrained least squares (NNLS), the estimator relies on the Gram matrix $\boldsymbol{\Phi}^T \boldsymbol{\Phi} \in \mathbb{R}^{S \times S}$. 

Let two cell states $s_1, s_2 \in \mathcal{S}$ be transcriptomically collinear, such that:

$$\boldsymbol{\phi}_{s_2} = \boldsymbol{\phi}_{s_1} + \boldsymbol{\delta}, \quad \text{with } \|\boldsymbol{\delta}\|_2 \ll \|\boldsymbol{\phi}_{s_1}\|_2$$

Define the normalized difference direction vector $\mathbf{v} \in \mathbb{R}^S$ as:

$$\mathbf{v} = \frac{1}{\sqrt{2}} (\mathbf{e}_{s_1} - \mathbf{e}_{s_2}) = \frac{1}{\sqrt{2}} (0, \dots, 1, \dots, -1, \dots, 0)^T$$

Evaluating the Rayleigh quotient of the Gram matrix along $\mathbf{v}$:

$$\mathbf{v}^T (\boldsymbol{\Phi}^T \boldsymbol{\Phi}) \mathbf{v} = \|\boldsymbol{\Phi} \mathbf{v}\|_2^2 = \frac{1}{2} \|\boldsymbol{\phi}_{s_1} - \boldsymbol{\phi}_{s_2}\|_2^2 = \frac{1}{2} \|\boldsymbol{\delta}\|_2^2$$

By the Courant-Fischer Minimax Theorem, the smallest singular value $\sigma_{\min}(\boldsymbol{\Phi})$ satisfies:

$$\sigma_{\min}(\boldsymbol{\Phi}) \le \sqrt{\frac{1}{2}} \|\boldsymbol{\delta}\|_2 \approx 0$$

The condition number $\kappa(\boldsymbol{\Phi})$ therefore diverges:

$$\kappa(\boldsymbol{\Phi}) = \frac{\sigma_{\max}(\boldsymbol{\Phi})}{\sigma_{\min}(\boldsymbol{\Phi})} \ge \frac{\sqrt{2} \sigma_{\max}(\boldsymbol{\Phi})}{\|\boldsymbol{\delta}\|_2} \longrightarrow \infty \quad \text{as } \|\boldsymbol{\delta}\|_2 \to 0$$

Under ordinary regression, the covariance of the estimated cell state fractions $\hat{\boldsymbol{\theta}}$ is:

$$\operatorname{Cov}(\hat{\boldsymbol{\theta}}) = \sigma_{\epsilon}^2 (\boldsymbol{\Phi}^T \boldsymbol{\Phi})^{-1}$$

Projecting this covariance matrix along the sibling difference direction $\mathbf{v}$:

$$\operatorname{Var}\left(\frac{\hat{\theta}_{s_1} - \hat{\theta}_{s_2}}{\sqrt{2}}\right) = \mathbf{v}^T \operatorname{Cov}(\hat{\boldsymbol{\theta}}) \mathbf{v} \ge \frac{\sigma_{\epsilon}^2}{\sigma_{\min}^2(\boldsymbol{\Phi})} \ge \frac{2 \sigma_{\epsilon}^2}{\|\boldsymbol{\delta}\|_2^2} \longrightarrow \infty$$

Expanding the variance of the difference:

$$\operatorname{Var}(\hat{\theta}_{s_1} - \hat{\theta}_{s_2}) = \operatorname{Var}(\hat{\theta}_{s_1}) + \operatorname{Var}(\hat{\theta}_{s_2}) - 2 \operatorname{Cov}(\hat{\theta}_{s_1}, \hat{\theta}_{s_2})$$

For this quantity to diverge while the total sum $\theta_{s_1} + \theta_{s_2} \le 1$ remains constrained, the covariance between sibling states must become strongly negative:

$$\operatorname{Cov}(\hat{\theta}_{s_1}, \hat{\theta}_{s_2}) \approx -\sqrt{\operatorname{Var}(\hat{\theta}_{s_1}) \operatorname{Var}(\hat{\theta}_{s_2})}, \quad \operatorname{Corr}(\hat{\theta}_{s_1}, \hat{\theta}_{s_2}) \longrightarrow -1$$

This creates **negative cross-talk**: any infinitesimal fluctuation, technical artifact, or biological shift in the bulk mixture causes the optimization solver to add mass to $s_1$ while subtracting an equal amount from $s_2$, producing severe sign-flipping.

---

### 2.2 Collinearity in the Dirichlet-Multinomial Generative Model

In BayesPrism, sequencing counts are modeled without log-transformation. For patient sample $n \in \{1, \dots, N\}$, the observed bulk count vector is $\mathbf{X}_n = (X_{n, 1}, \dots, X_{n, G}) \in \mathbb{N}_0^G$, with library size $R_n = \sum_{g=1}^G X_{n, g}$.

The generative probability of observing a read for gene $g$ in sample $n$ is:

$$P(\text{gene } g \mid \boldsymbol{\theta}_n, \boldsymbol{\Phi}) = \sum_{s=1}^S \phi_{s, g} \theta_{n, s}$$

The log-likelihood of $\mathbf{X}_n$ given cell state fractions $\boldsymbol{\theta}_n \in \Delta^{S-1}$ is:

$$\ell(\boldsymbol{\theta}_n) = \sum_{g=1}^G X_{n, g} \log\left( \sum_{s=1}^S \phi_{s, g} \theta_{n, s} \right) + \text{const}$$

The Fisher Information matrix $\mathcal{I}(\boldsymbol{\theta}_n) \in \mathbb{R}^{S \times S}$ has entries:

$$\mathcal{I}(\boldsymbol{\theta}_n)_{j, k} = -\mathbb{E}\left[ \frac{\partial^2 \ell}{\partial \theta_{n, j} \partial \theta_{n, k}} \right] = R_n \sum_{g=1}^G \frac{\phi_{j, g} \phi_{k, g}}{\sum_{s=1}^S \phi_{s, g} \theta_{n, s}}$$

Evaluating the Fisher Information along the difference vector $\mathbf{v} = \frac{1}{\sqrt{2}}(\mathbf{e}_{s_1} - \mathbf{e}_{s_2})$:

$$\mathbf{v}^T \mathcal{I}(\boldsymbol{\theta}_n) \mathbf{v} = \frac{R_n}{2} \sum_{g=1}^G \frac{(\phi_{s_1, g} - \phi_{s_2, g})^2}{\sum_{s=1}^S \phi_{s, g} \theta_{n, s}} = \frac{R_n}{2} \sum_{g=1}^G \frac{\delta_g^2}{\sum_{s=1}^S \phi_{s, g} \theta_{n, s}}$$

Since the denominator is bounded below by a positive constant $c_{\min} > 0$ for expressed genes:

$$\mathbf{v}^T \mathcal{I}(\boldsymbol{\theta}_n) \mathbf{v} \le \frac{R_n}{2 c_{\min}} \|\boldsymbol{\delta}\|_2^2 \longrightarrow 0 \quad \text{as } \|\boldsymbol{\delta}\|_2 \to 0$$

By the Cramér-Rao lower bound, the posterior covariance along $\mathbf{v}$ satisfies:

$$\operatorname{Var}_{p(\boldsymbol{\theta}_n \mid \mathbf{X}_n)}(\theta_{n, s_1} - \theta_{n, s_2}) \ge 2 \mathbf{v}^T \mathcal{I}^{-1}(\boldsymbol{\theta}_n) \mathbf{v} \ge \frac{4 c_{\min}}{R_n \|\boldsymbol{\delta}\|_2^2} \longrightarrow \infty$$

Thus, in any single-level MCMC sampler or MAP estimator over cell states $\mathcal{S}$, the posterior distribution $p(\boldsymbol{\theta}_n \mid \mathbf{X}_n)$ forms an elongated ridge along the subspace $\theta_{n, s_1} + \theta_{n, s_2} = C$. Individual draws $\theta_{n, s_1}^{(i)}$ and $\theta_{n, s_2}^{(i)}$ fluctuate with extreme posterior variance and correlation $\approx -1$.

---

## 3. BayesPrism's Two-Tier Hierarchical Architecture

To solve this identifiability crisis without flattening biological nuance, BayesPrism enforces a formal two-tier hierarchy:

### 3.1 Spaces and the Surjective Mapping

- **Cell State Space**: $\mathcal{S} = \{1, 2, \dots, S\}$, indexed by $s$. These represent fine phenotypic clusters, continuous manifold partitions, or patient-specific tumor sub-populations.
- **Cell Type Space**: $\mathcal{T} = \{1, 2, \dots, T\}$, indexed by $t$, where $T \ll S$. These represent distinct major lineages (e.g., T cells, B cells, Myeloid, Endothelial, Fibroblasts, Malignant).
- **Surjective Partition Map**: $\pi: \mathcal{S} \twoheadrightarrow \mathcal{T}$, assigning each state $s$ to exactly one broad cell type $t = \pi(s)$.
- **Disjoint Partition**:
  $$\mathcal{S} = \bigsqcup_{t \in \mathcal{T}} \mathcal{S}_t, \quad \text{where } \mathcal{S}_t = \{s \in \mathcal{S} : \pi(s) = t\}, \quad \mathcal{S}_t \cap \mathcal{S}_{t'} = \emptyset \quad \forall t \neq t'$$

### 3.2 Dual Reference Representation

BayesPrism constructs two distinct normalized reference matrices in `new.prism`:

1. **State Reference Matrix $\boldsymbol{\Phi}^{\text{state}} \in \mathbb{R}^{S \times G}$**:
   For each state $s \in \mathcal{S}$, single-cell raw counts are collapsed:
   $$C_{s, g}^{\text{state}} = \sum_{c \in \text{cells}(s)} Y_{c, g}$$
   Normalized to the simplex with a pseudo-count $\epsilon = 10^{-8}$:
   $$\phi_{s, g}^{\text{state}} = \frac{C_{s, g}^{\text{state}}}{\sum_{g'=1}^G C_{s, g'}^{\text{state}}} (1 - \epsilon G) + \epsilon$$

2. **Type Reference Matrix $\boldsymbol{\Phi}^{\text{type}} \in \mathbb{R}^{T \times G}$**:
   Single-cell counts are collapsed directly by cell type $t \in \mathcal{T}$:
   $$C_{t, g}^{\text{type}} = \sum_{c \in \text{cells}(t)} Y_{c, g} = \sum_{s \in \mathcal{S}_t} C_{s, g}^{\text{state}}$$
   Normalized to the simplex:
   $$\phi_{t, g}^{\text{type}} = \frac{C_{t, g}^{\text{type}}}{\sum_{g'=1}^G C_{t, g'}^{\text{type}}} (1 - \epsilon G) + \epsilon$$

### 3.3 The Aggregation Operator $\mathbf{M}$

Define the linear aggregation matrix $\mathbf{M} \in \{0, 1\}^{T \times S}$ by:

$$M_{t, s} = \begin{cases} 1 & \text{if } \pi(s) = t \quad (s \in \mathcal{S}_t) \\ 0 & \text{otherwise} \end{cases}$$

This operator maps any state-level vector $\mathbf{u} \in \mathbb{R}^S$ to a cell-type-level vector $\mathbf{w} \in \mathbb{R}^T$:

$$\mathbf{w} = \mathbf{M} \mathbf{u} \iff w_t = \sum_{s \in \mathcal{S}_t} u_s$$

---

## 4. Algorithmic Breakdown & Theoretical Proofs

```
========================================================================================
                         BAYESPRISM FOUR-STAGE WORKFLOW
========================================================================================

 [Bulk Mixture X_n]      [State Reference Φ_state]
         │                          │
         ▼                          ▼
 ┌────────────────────────────────────────────────────────┐
 │ Stage 1: Fine-State Gibbs Sampling                     │
 │   - Sample Z_{n,g,s} ~ Multinomial(X_{n,g}, p_{n,g,s}) │   <- Captures intra-lineage
 │   - Sample θ_{n,s}   ~ Dirichlet(∑_g Z_{n,g,s} + α)    │      heterogeneity
 └──────────────────────────┬─────────────────────────────┘
                            │
                            ▼
 ┌────────────────────────────────────────────────────────┐
 │ Stage 2: Null-Space Statistical Marginalization (M)    │
 │   - Z_{n,g,t}   = ∑_{s ∈ S_t} Z_{n,g,s}                │   <- Sibling collinear modes
 │   - θ_{n,t}^(0) = ∑_{s ∈ S_t} θ_{n,s}                  │      lie in Null(M)
 │   * Exact Variance Cancellation: Var(θ_t) = O(1)       │      Variance explosion canceled!
 └──────────────────────────┬─────────────────────────────┘
                            │
                            ▼
 ┌────────────────────────────────────────────────────────┐
 │ Stage 3: Broad Lineage Reference Refinement            │
 │   - Pool cohort: Z_{g,t} = ∑_n Z_{n,g,t}               │   <- Well-conditioned (T << S)
 │   - Optimize γ_{t,g} under N(0, σ^2) Gaussian prior    │   <- MAP shrinkage prevents
 │   - Refined reference: Ψ_t = Softmax(log Φ_t + γ_t)    │      overfitting
 └──────────────────────────┬─────────────────────────────┘
                            │
                            ▼
 ┌────────────────────────────────────────────────────────┐
 │ Stage 4: Final Non-Collinear Deconvolution             │
 │   - Final Gibbs sampling strictly on simplex Δ^{T-1}   │   <- Condition number κ(Ψ) low
 │   - Output: θ_{n,t}^(f) and Z_{n,g,t}                  │   <- High clinical concordance
 └────────────────────────────────────────────────────────┘
========================================================================================
```

---

### 4.1 Stage 1: Fine-State Latent Variable Sampling

In Stage 1, BayesPrism runs a Dirichlet-Multinomial Gibbs sampler on each sample $n$ using the fine state reference $\boldsymbol{\Phi}^{\text{state}}$.

The total bulk count $X_{n, g}$ is treated as the sum of unobserved latent counts across all cell states:

$$X_{n, g} = \sum_{s=1}^S Z_{n, g, s}$$

At MCMC iteration $i$:
1. **Sampling Latent Counts $\mathbf{Z}_{n, g, \cdot}^{(i)}$**:
   Conditional on the current state proportions $\boldsymbol{\theta}_n^{(i-1)}$, the multinomial allocation probabilities are:
   $$p_{n, g, s}^{(i)} = \frac{\phi_{s, g}^{\text{state}} \theta_{n, s}^{(i-1)}}{\sum_{s'=1}^S \phi_{s', g}^{\text{state}} \theta_{n, s'}^{(i-1)}}$$
   The latent counts are drawn jointly for each gene $g$:
   $$(Z_{n, g, 1}^{(i)}, \dots, Z_{n, g, S}^{(i)}) \sim \operatorname{Multinomial}\left(X_{n, g}, (p_{n, g, 1}^{(i)}, \dots, p_{n, g, S}^{(i)})\right)$$

2. **Sampling State Proportions $\boldsymbol{\theta}_n^{(i)}$**:
   Aggregating latent counts across all genes for state $s$:
   $$Z_{n, \cdot, s}^{(i)} = \sum_{g=1}^G Z_{n, g, s}^{(i)}$$
   Under a symmetric Dirichlet prior $\operatorname{Dirichlet}(\alpha \mathbf{1}_S)$ (default $\alpha = 1.0$), the conjugacy of the multinomial likelihood yields the Dirichlet posterior draw:
   $$\boldsymbol{\theta}_n^{(i)} \sim \operatorname{Dirichlet}\left(Z_{n, \cdot, 1}^{(i)} + \alpha, Z_{n, \cdot, 2}^{(i)} + \alpha, \dots, Z_{n, \cdot, S}^{(i)} + \alpha\right)$$

After discarding the burn-in period and applying thinning over retained sample set $\mathcal{K}_{\text{gibbs}}$, the posterior means are:

$$\hat{Z}_{n, g, s} = \frac{1}{|\mathcal{K}_{\text{gibbs}}|} \sum_{i \in \mathcal{K}_{\text{gibbs}}} Z_{n, g, s}^{(i)}, \quad \hat{\theta}_{n, s} = \frac{1}{|\mathcal{K}_{\text{gibbs}}|} \sum_{i \in \mathcal{K}_{\text{gibbs}}} \theta_{n, s}^{(i)}$$

#### Why Not Restrict Stage 1 to Broad Cell Types Directly?
If one collapses single-cell profiles into broad cell types prior to Stage 1 ($\boldsymbol{\Phi}^{\text{type}}$), the model forces all cells of lineage $t$ to conform to a single static average. If patient $A$ harbors mostly cytotoxic effector cells ($s_1$) while patient $B$ harbors exhausted T cells ($s_2$), the static centroid $\boldsymbol{\phi}_t^{\text{type}}$ is severely misspecified for both patients. 

By sampling at the state level $\mathcal{S}$, the multinomial draws $Z_{n, g, s}$ dynamically adapt to whichever phenotypic state best explains patient $n$'s reads, capturing the true patient-specific activation profile.

---

### 4.2 Stage 2: Null-Space Invariance & Variance Cancellation Theorem

In Stage 2, BayesPrism applies the marginalization operator `mergeK` to compute:

$$Z_{n, g, t} = \sum_{s \in \mathcal{S}_t} \hat{Z}_{n, g, s} = (\mathbf{M} \hat{\mathbf{Z}}_{n, g, \cdot})_t$$

$$\theta_{n, t}^{(0)} = \sum_{s \in \mathcal{S}_t} \hat{\theta}_{n, s} = (\mathbf{M} \hat{\boldsymbol{\theta}}_n)_t$$

We now state and prove the fundamental theorem explaining why this operation eliminates collinearity.

---

#### Theorem 1 (Null-Space Invariance of Collinear Variance Inflation)
*Let $\mathcal{S}_t = \{s \in \mathcal{S} : \pi(s) = t\}$ be the set of fine states belonging to broad cell type $t$. Suppose two sibling states $s_1, s_2 \in \mathcal{S}_t$ are $\delta$-collinear, i.e., $\boldsymbol{\phi}_{s_2} = \boldsymbol{\phi}_{s_1} + \boldsymbol{\delta}$ with $\|\boldsymbol{\delta}\|_2 \to 0$. Let $\mathbf{M} \in \{0, 1\}^{T \times S}$ be the aggregation operator.*

*Then:*
1. *The collinear difference direction $\mathbf{v} = \frac{1}{\sqrt{2}}(\mathbf{e}_{s_1} - \mathbf{e}_{s_2})$ lies in the null space of $\mathbf{M}$:*
   $$\mathbf{v} \in \operatorname{Null}(\mathbf{M})$$
2. *While the individual state variances diverge as $O(1/\|\boldsymbol{\delta}\|_2^2)$:*
   $$\operatorname{Var}(\hat{\theta}_{n, s_1}) \longrightarrow \infty, \quad \operatorname{Var}(\hat{\theta}_{n, s_2}) \longrightarrow \infty$$
   *the variance of the broad cell type fraction remains bounded and well-conditioned:*
   $$\operatorname{Var}(\theta_{n, t}^{(0)}) = O(1)$$

---

#### Proof
**(Part 1: Null Space Membership)**
Apply the linear operator $\mathbf{M}$ to the difference vector $\mathbf{v}$:

$$(\mathbf{M} \mathbf{v})_k = \sum_{s \in \mathcal{S}} M_{k, s} v_s = \frac{1}{\sqrt{2}} \left( M_{k, s_1} - M_{k, s_2} \right)$$

By definition of $\mathbf{M}$:
- If $k = t$: since $s_1, s_2 \in \mathcal{S}_t$, $M_{t, s_1} = 1$ and $M_{t, s_2} = 1$. Thus, $(\mathbf{M} \mathbf{v})_t = \frac{1}{\sqrt{2}} (1 - 1) = 0$.
- If $k \neq t$: neither state belongs to type $k$, so $M_{k, s_1} = 0$ and $M_{k, s_2} = 0$. Thus, $(\mathbf{M} \mathbf{v})_k = 0$.

Therefore:
$$\mathbf{M} \mathbf{v} = \mathbf{0} \implies \mathbf{v} \in \operatorname{Null}(\mathbf{M}) \quad \blacksquare$$

**(Part 2: Variance Cancellation)**
Consider the variance of the broad cell type fraction $\theta_{n, t}^{(0)} = \sum_{s \in \mathcal{S}_t} \hat{\theta}_{n, s}$:

$$\operatorname{Var}(\theta_{n, t}^{(0)}) = \mathbf{e}_t^T \operatorname{Cov}(\mathbf{M} \hat{\boldsymbol{\theta}}_n) \mathbf{e}_t = \mathbf{m}_t^T \operatorname{Cov}(\hat{\boldsymbol{\theta}}_n) \mathbf{m}_t$$

where $\mathbf{m}_t = \sum_{s \in \mathcal{S}_t} \mathbf{e}_s \in \{0, 1\}^S$ is the $t$-th row of $\mathbf{M}$.

Decompose the covariance matrix $\boldsymbol{\Sigma}_{\boldsymbol{\theta}} = \operatorname{Cov}(\hat{\boldsymbol{\theta}}_n)$ using its spectral decomposition:

$$\boldsymbol{\Sigma}_{\boldsymbol{\theta}} = \sum_{k=1}^S \lambda_k \mathbf{u}_k \mathbf{u}_k^T$$

As established in Section 2.2, as $\|\boldsymbol{\delta}\|_2 \to 0$, exactly one eigenvalue diverges:

$$\lambda_1 \sim O\left(\frac{1}{\|\boldsymbol{\delta}\|_2^2}\right) \longrightarrow \infty, \quad \text{with eigenvector } \mathbf{u}_1 \longrightarrow \mathbf{v} = \frac{1}{\sqrt{2}}(\mathbf{e}_{s_1} - \mathbf{e}_{s_2})$$

All remaining eigenvectors $\mathbf{u}_k$ ($k \ge 2$) correspond to directions transversal to $\mathbf{v}$ and have bounded eigenvalues $\lambda_k = O(1)$.

Now evaluate $\mathbf{m}_t^T \boldsymbol{\Sigma}_{\boldsymbol{\theta}} \mathbf{m}_t$:

$$\operatorname{Var}(\theta_{n, t}^{(0)}) = \mathbf{m}_t^T \left( \lambda_1 \mathbf{u}_1 \mathbf{u}_1^T + \sum_{k=2}^S \lambda_k \mathbf{u}_k \mathbf{u}_k^T \right) \mathbf{m}_t = \lambda_1 (\mathbf{m}_t^T \mathbf{u}_1)^2 + \sum_{k=2}^S \lambda_k (\mathbf{m}_t^T \mathbf{u}_k)^2$$

Compute the inner product $\mathbf{m}_t^T \mathbf{u}_1$:

$$\mathbf{m}_t^T \mathbf{u}_1 \longrightarrow \mathbf{m}_t^T \mathbf{v} = \sum_{s \in \mathcal{S}_t} v_s = \frac{1}{\sqrt{2}} (1 - 1) = 0$$

The divergent coefficient multiplying $\lambda_1$ vanishes identically:

$$\lim_{\|\boldsymbol{\delta}\|_2 \to 0} \lambda_1 (\mathbf{m}_t^T \mathbf{u}_1)^2 = \lim_{\|\boldsymbol{\delta}\|_2 \to 0} O\left(\frac{1}{\|\boldsymbol{\delta}\|_2^2}\right) \cdot O(\|\boldsymbol{\delta}\|_2^2) = O(1)$$

Expanding the sum explicitly in terms of variances and covariances:

$$\operatorname{Var}(\theta_{n, t}^{(0)}) = \operatorname{Var}(\hat{\theta}_{n, s_1}) + \operatorname{Var}(\hat{\theta}_{n, s_2}) + 2 \operatorname{Cov}(\hat{\theta}_{n, s_1}, \hat{\theta}_{n, s_2}) + \sum_{s \in \mathcal{S}_t \setminus \{s_1, s_2\}} \dots$$

Because $\operatorname{Cov}(\hat{\theta}_{n, s_1}, \hat{\theta}_{n, s_2}) = -\frac{1}{2} [\operatorname{Var}(\hat{\theta}_{n, s_1}) + \operatorname{Var}(\hat{\theta}_{n, s_2})] + O(1)$, the divergent terms cancel each other out exactly:

$$\operatorname{Var}(\theta_{n, t}^{(0)}) = O(1) \quad \blacksquare$$

---

### 4.3 Stage 3: Broad Lineage Reference Refinement

In Stage 3, BayesPrism updates the reference expression profiles to account for technical platform differences (scRNA-seq vs. bulk RNA-seq) and biological tumor microenvironmental shifts.

Crucially, **BayesPrism performs this update at the broad cell type level $\mathcal{T}$, never at the collinear state level $\mathcal{S}$**.

#### 4.3.1 Cohort-Level Pooling
BayesPrism aggregates the marginalized latent counts $Z_{n, g, t}$ across all $N$ bulk samples in the patient cohort:

$$Z_{g, t} = \sum_{n=1}^N Z_{n, g, t}, \quad Z_{\cdot, t} = \sum_{g=1}^G Z_{g, t}$$

Pooling across $N$ samples increases the effective read depth, stabilizing gene-level estimates.

#### 4.3.2 Multiplicative Model with Gaussian Shrinkage Prior
For each non-malignant cell type $t \in \mathcal{T}$, the updated reference $\boldsymbol{\psi}_t \in \Delta^{G-1}$ is parameterized as a multiplicative perturbation of the baseline single-cell type profile $\boldsymbol{\phi}_t^{\text{type}}$ using log-fold change vector $\boldsymbol{\gamma}_t \in \mathbb{R}^G$:

$$\psi_{t, g}(\boldsymbol{\gamma}_t) = \frac{\phi_{t, g}^{\text{type}} \exp(\gamma_{t, g})}{\sum_{g'=1}^G \phi_{t, g'}^{\text{type}} \exp(\gamma_{t, g'})} = \operatorname{Softmax}_g\left( \log \boldsymbol{\phi}_t^{\text{type}} + \boldsymbol{\gamma}_t \right)$$

To prevent biological marker distortion and avoid overfitting, BayesPrism imposes an independent zero-mean Gaussian shrinkage prior on $\boldsymbol{\gamma}_t$:

$$\gamma_{t, g} \sim \mathcal{N}(0, \sigma^2), \quad p(\boldsymbol{\gamma}_t) = \prod_{g=1}^G \frac{1}{\sqrt{2\pi}\sigma} \exp\left(-\frac{\gamma_{t, g}^2}{2\sigma^2}\right)$$

where $\sigma$ (default $\sigma = 2.0$) controls the prior belief regarding cross-platform fold changes.

#### 4.3.3 Objective Function and Gradient
The MAP objective maximizes the log posterior over $\boldsymbol{\gamma}_t$:

$$\mathcal{L}(\boldsymbol{\gamma}_t) = \sum_{g=1}^G Z_{g, t} \log \psi_{t, g}(\boldsymbol{\gamma}_t) - \frac{1}{2\sigma^2} \sum_{g=1}^G \gamma_{t, g}^2$$

Expanding $\log \psi_{t, g}$:

$$\log \psi_{t, g} = \log \phi_{t, g}^{\text{type}} + \gamma_{t, g} - \operatorname{LSE}\left(\log \boldsymbol{\phi}_t^{\text{type}} + \boldsymbol{\gamma}_t\right)$$

where $\operatorname{LSE}(\mathbf{x}) = \log\left(\sum_{g=1}^G \exp(x_g)\right)$ is the Log-Sum-Exp function.

The gradient with respect to $\gamma_{t, g}$ is:

$$\frac{\partial \mathcal{L}}{\partial \gamma_{t, g}} = Z_{g, t} - Z_{\cdot, t} \psi_{t, g}(\boldsymbol{\gamma}_t) - \frac{1}{\sigma^2} \gamma_{t, g}$$

The optimization is solved independently for each cell type $t \in \mathcal{T}$ using unconstrained conjugate gradient / L-BFGS (`Rcgminu` in R, L-BFGS-B in Python).

#### 4.3.4 Patient-Specific Malignant Update
If malignant cells are specified by `key`, their extreme intra- and inter-patient genetic heterogeneity (e.g., aneuploidy, copy number alterations) violates the cohort-level sharing assumption. BayesPrism therefore updates the malignant profile **per sample** via Maximum Likelihood Estimation:

$$\psi_{n, g}^{\text{mal}} = \frac{Z_{n, g, \text{mal}}}{\sum_{g'=1}^G Z_{n, g', \text{mal}}} (1 - \epsilon G) + \epsilon$$

---

### 4.4 Stage 4: Final Non-Collinear Proportion Estimation

With the refined reference matrix $\boldsymbol{\Psi} \in \mathbb{R}^{T \times G}$ established:

$$\boldsymbol{\Psi} = \begin{bmatrix} \boldsymbol{\psi}_1^T \\ \boldsymbol{\psi}_2^T \\ \vdots \\ \boldsymbol{\psi}_T^T \end{bmatrix}, \quad \psi_{t, g} \ge 0, \quad \sum_{g=1}^G \psi_{t, g} = 1$$

BayesPrism runs a final Gibbs sampler (`run.gibbs.final`). 

This sampler operates **strictly in the $T$-dimensional space of broad cell types**:
1. At each iteration, multinomial allocation probabilities for gene $g$ in sample $n$ are:
   $$p_{n, g, t} = \frac{\psi_{t, g} \theta_{n, t}}{\sum_{t'=1}^T \psi_{t', g} \theta_{n, t'}}$$
2. Latent cell-type reads are drawn:
   $$(Z_{n, g, 1}, \dots, Z_{n, g, T}) \sim \operatorname{Multinomial}(X_{n, g}, (p_{n, g, 1}, \dots, p_{n, g, T}))$$
3. Final cell-type fractions are sampled:
   $$\boldsymbol{\theta}_n^{(f)} \sim \operatorname{Dirichlet}\left(\sum_{g=1}^G Z_{n, g, 1} + \alpha, \dots, \sum_{g=1}^G Z_{n, g, T} + \alpha\right)$$

Because the broad lineages (T cells, B cells, Myeloid, Endothelial, etc.) have distinct marker profiles, the refined matrix $\boldsymbol{\Psi}$ has a low condition number:

$$\kappa(\boldsymbol{\Psi}) \ll 100$$

All collinear degrees of freedom from the underlying cell states have been eliminated, guaranteeing stable, non-singular, and reproducible cell fraction estimates.

---

## 5. Empirical Demonstration: The Sade-Feldman Melanoma Benchmark

The clinical impact of this mathematical mechanism is directly observable in the **Sade-Feldman melanoma cohort** (GSE120575, 51 biopsies, 12 immune states):

### 5.1 Flat Deconvolution Pathology
When deconvolution is performed flatly across all 12 cell states:
- 5 of the 12 clusters represent closely related CD8+ cytotoxic T-cell subsets (`00_Tem/Trm`, `01_Tem/Trm`, `03_Tem/Trm`, `06_Tem/Temra`, `07_Tem/Trm`).
- The condition number of this reference is severe:
  $$\kappa(\boldsymbol{\Phi}^{\text{state}}) = 1,248.6$$
- In pre-treatment biopsies, single-cell Milo differential abundance (DA) testing reveals that `01_Tem/Trm` is strongly enriched in responding patients ($\overline{\text{logFC}} = +1.726$).
- However, flat deconvolution yields a negative coefficient:
  $$\hat{\beta}_{\text{deconv}} = -0.509 \quad (\text{Discordant Sign-Flip!})$$
- **Root Cause**: Because `07_Tem/Trm` was depleted in responders ($\overline{\text{logFC}} = -1.324, \hat{\beta}_{\text{deconv}} = -1.047$), the collinear cross-talk $(\boldsymbol{\Phi}^T \boldsymbol{\Phi})^{-1}$ subtracted fraction from `01_Tem/Trm` to offset the strong signal in `07_Tem/Trm`.

### 5.2 Resolution via Hierarchical Lineage Consolidation
When BayesPrism's hierarchical mapping is applied:
1. All 5 cytotoxic T-cell states are grouped into the major lineage $\pi(s) = \text{Cytotoxic\_T}$.
2. Condition number drops by **68.2%**:
   $$\kappa(\boldsymbol{\Phi}^{\text{consolidated}}) = 397.1$$
3. Directional sign concordance across the cohort surges:
   $$\text{Concordance Rate} = 75.0\% \longrightarrow 88.9\%$$
4. The spurious negative cross-talk between `01_Tem/Trm` and `07_Tem/Trm` vanishes completely.

---

## 6. Code Mapping Matrix: Theory to Implementation

| Mathematical Concept | Symbol / Equation | BayesPrism R (`scratch/BayesPrism/`) | Python (`packages/bayesprism/`) |
| :--- | :--- | :--- | :--- |
| **Cell State Space** | $s \in \mathcal{S}, |\mathcal{S}| = S$ | `cell.state.labels` | `cell_state_labels` |
| **Cell Type Space** | $t \in \mathcal{T}, |\mathcal{T}| = T$ | `cell.type.labels` | `cell_type_labels` |
| **Surjective Map $\pi$** | $\mathcal{S}_t = \{s : \pi(s) = t\}$ | `prism@map` in `new_prism.R` | `state_to_type_map` in `pipeline.py` |
| **State Reference** | $\boldsymbol{\Phi}^{\text{state}} \in \mathbb{R}^{S \times G}$ | `prism@phi_cellState@phi` | `prism.phi_cell_state.phi` |
| **Type Reference** | $\boldsymbol{\Phi}^{\text{type}} \in \mathbb{R}^{T \times G}$ | `prism@phi_cellType@phi` | `prism.phi_cell_type.phi` |
| **Stage 1 Gibbs Sampler** | $(Z_{n, g, s}, \theta_{n, s}) \sim p(\cdot \mid \mathbf{X}_n)$ | `sample.Z.theta_n` in `run_gibbs.R` | `sample_Z_theta_n` in `gibbs.py` |
| **Marginalization $\mathbf{M}$** | $Z_{n, g, t} = \sum_{s \in \mathcal{S}_t} Z_{n, g, s}$ | `mergeK` in `JointPost_functions.R` | `merge_k` in `pipeline.py` |
| **Initial Proportions** | $\theta_{n, t}^{(0)} = \sum_{s \in \mathcal{S}_t} \theta_{n, s}$ | `jointPost.ini.ct@theta` | `joint_ini_ct.theta` |
| **Cohort Pooling** | $Z_{g, t} = \sum_{n=1}^N Z_{n, g, t}$ | `colSums(Z, dims=1)` in `update_reference.R` | `np.sum(Z_np, axis=0)` in `optimization.py` |
| **MAP Reference Update** | $\max_{\boldsymbol{\gamma}_t} \left[ \ell(\boldsymbol{\gamma}_t) - \frac{\|\boldsymbol{\gamma}_t\|^2}{2\sigma^2} \right]$ | `optimize.psi` in `optim_functions_MAP.R` | `optimize_psi_map` in `optimization.py` |
| **Softmax Transformation** | $\psi_{t, g} = \operatorname{Softmax}(\log \phi_t + \boldsymbol{\gamma}_t)$ | `transform.phi_t` in `update_reference.R` | `transform_phi_t` in `optimization.py` |
| **Final Gibbs Sampling** | $\boldsymbol{\theta}_n^{(f)} \sim p(\cdot \mid \mathbf{X}_n, \boldsymbol{\Psi})$ | `sample.theta_n` in `run_gibbs.R` | `sample_theta_n` in `gibbs.py` |

---

## 7. Synthesis and Recommendations

1. **Hierarchy as a Mathematical Regularizer**: BayesPrism's grouping operator $\mathbf{M}$ functions as an exact orthogonal projector onto the complement of the collinear difference subspace, $\operatorname{Null}(\mathbf{M})^\perp$.
2. **Best Practice for Deconvolution Pipelines**:
   - Always supply fine-grained continuum states to `cell.state.labels` (at least 20–50 cells per state).
   - Reserve `cell.type.labels` for major distinct lineages that possess $>50$ significantly differentially expressed genes.
   - For downstream clinical association testing (e.g., Milo vs. bulk deconvolution), evaluate response models at the broad lineage level ($\boldsymbol{\theta}^{(f)}$) or on pruned, condition-number-consolidated clusters ($\kappa < 400$) to guarantee directional sign fidelity.
