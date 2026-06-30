# Mathematical Explanation of `test_difference_in_means` in Selective Inference

This document explains the mathematical theory and implementation details behind the `test_difference_in_means` functions used for post-selection inference in clustering.

---

## 1. The Core Problem: Selection Bias ("Double Dipping")

When we apply a clustering algorithm (such as a Gaussian Mixture Model or K-Means) to a dataset $X \in \mathbb{R}^{n \times d}$ and then test whether two clusters $c_1$ and $c_2$ have different means along the direction of their observed difference, we are reusing the same data for two distinct purposes:
1. **Selection:** Defining the clusters $c_1, c_2$ and the test direction.
2. **Inference:** Testing the statistical significance of the difference in means.

Under the null hypothesis $H_0$ that the true cluster means are equal, a naive hypothesis test (e.g., a standard $t$-test or $z$-test) yields highly biased $p$-values that are shifted toward $0$. This occurs because the clustering algorithm groups similar points together, artificially inflating the difference in sample means even when the underlying data is purely random noise.

To obtain valid $p$-values, we must perform **Selective Inference** (post-selection inference) by conditioning the test statistic on the event that the clustering algorithm would choose the exact same partition of the data.

---

## 2. Mathematical Setup

Let $X \in \mathbb{R}^{n \times d}$ be the observed data matrix where each row $X_i \in \mathbb{R}^d$ is an independent observation, modeled as:
$$X \sim \mathcal{N}(M, \sigma^2 I_{n \times d})$$
where $M \in \mathbb{R}^{n \times d}$ is the matrix of true means and $\sigma^2$ is the noise variance.

Let $C(X) = \hat{y} \in \{1, \dots, K\}^n$ be the cluster labels assigned to the data points. Let $C_1, C_2 \subset \{1, \dots, n\}$ be the index sets of the two clusters being compared, with sizes $n_1 = |C_1|$ and $n_2 = |C_2|$.

### 2.1 The Direction of Interest
The sample means of the two clusters are:
$$\bar{X}_1 = \frac{1}{n_1} \sum_{i \in C_1} X_i, \quad \bar{X}_2 = \frac{1}{n_2} \sum_{i \in C_2} X_i$$

We define the unit difference vector along the direction of the observed difference:
$$\delta = \frac{\bar{X}_1 - \bar{X}_2}{\|\bar{X}_1 - \bar{X}_2\|_2} \in \mathbb{R}^d$$

### 2.2 The Contrast Matrix
We construct a projection matrix $V \in \mathbb{R}^{n \times d}$ that isolates the difference between the two sample means:
$$V_i = \begin{cases} \frac{\delta}{n_1} & \text{if } i \in C_1 \\ -\frac{\delta}{n_2} & \text{if } i \in C_2 \\ 0 & \text{otherwise} \end{cases}$$

### 2.3 The Observed Test Statistic
The observed test statistic $\phi_{\text{obs}}$ is the projection of $X$ along $V$:
$$\phi(X) = \langle V, X \rangle_F = \text{tr}(V^T X) = \sum_{i \in C_1} \frac{\delta^T X_i}{n_1} - \sum_{i \in C_2} \frac{\delta^T X_i}{n_2} = \delta^T (\bar{X}_1 - \bar{X}_2) = \|\bar{X}_1 - \bar{X}_2\|_2$$

Under the null hypothesis $H_0: \text{tr}(V^T M) = 0$, the statistic $\phi(X) \sim \mathcal{N}(0, \sigma^2 \|V\|_F^2)$ without conditioning, where:
$$\text{Var}(\phi(X)) = \sigma^2 \|V\|_F^2 = \sigma^2 \left( \frac{1}{n_1} + \frac{1}{n_2} \right)$$
which yields the naive standard deviation $\sigma_{\phi} = \sigma \sqrt{\|V\|_F^2}$.

---

## 3. Conditioning and 1D Perturbation

To eliminate the selection bias, we condition on the selection event $C(X) = \hat{y}_{\text{obs}}$, as well as the components of the data orthogonal to the direction of interest:
$$X_{\text{perp}} = P_{V^\perp} X = X - \frac{\text{tr}(V^T X)}{\|V\|_F^2} V$$

Conditioning on $P_{V^\perp} X$ restricts the stochasticity to a 1D line in the space of datasets, parameterized by $z \in \mathbb{R}$:
$$X(z) = X + (z - \phi_{\text{obs}}) V$$
Note that at $z = \phi_{\text{obs}}$, we recover the original dataset $X(\phi_{\text{obs}}) = X$. 

Under $H_0$, the conditional distribution of $z$ given the selection event is:
$$z \mid \{C(X(z)) = \hat{y}_{\text{obs}}\} \;\sim\; \mathcal{N}\left(0, \sigma^2 \|V\|_F^2\right) \text{ restricted to the support } \mathcal{S} = \{z \in \mathbb{R} \mid C(X(z)) = \hat{y}_{\text{obs}}\}$$

---

## 4. Implementation 1: GMM MCMC Selective Inference

For a Gaussian Mixture Model (GMM), the partition function $C(X(z))$ is determined via the Expectation-Maximization (EM) algorithm, which is highly non-linear and sensitive to initialization. Consequently, the set $\mathcal{S}$ has no analytical form.

In [gmm_mcmc_inference.py](file:///Users/halu/Code/tme_analysis/scripts/selective_inference/gmm_mcmc_inference.py#L20), the conditional distribution is sampled using a **Metropolis-Hastings Hit-and-Run sampler**:

1. **Initial State:** Start at the observed test statistic: $z_{\text{current}} = \phi_{\text{obs}}$.
2. **Proposal:** Propose a new candidate $z_{\text{prop}} \sim \mathcal{N}(z_{\text{current}}, \tau^2)$, where the step size $\tau$ is tuned to the null standard deviation ($\tau = 0.5 \sigma_{\phi}$).
3. **Perturbation and Clustering:** Construct the perturbed dataset:
   $$X_{\text{prop}} = X + (z_{\text{prop}} - \phi_{\text{obs}}) V$$
   Fit a new GMM (with a fixed random state to ensure deterministic labels) and obtain $\hat{y}_{\text{prop}}$.
4. **Metropolis-Hastings Decision:**
   - If the clustering matches the original labels ($\hat{y}_{\text{prop}} = \hat{y}_{\text{obs}}$), calculate the acceptance probability $\alpha$:
     $$\alpha = \min\left(1, \exp\left(-\frac{1}{2 \sigma^2 \|V\|_F^2} (z_{\text{prop}}^2 - z_{\text{current}}^2)\right)\right)$$
     In log space, this is:
     $$\log \alpha = -0.5 \left( \left(\frac{z_{\text{prop}}}{\sigma_{\phi}}\right)^2 - \left(\frac{z_{\text{current}}}{\sigma_{\phi}}\right)^2 \right)$$
     With probability $\alpha$, update $z_{\text{current}} = z_{\text{prop}}$.
   - If $\hat{y}_{\text{prop}} \neq \hat{y}_{\text{obs}}$, the proposal falls outside the support $\mathcal{S}$, and we reject it (acceptance probability 0).
5. **Selective $p$-value:** After a burn-in period, the empirical selective $p$-value is computed as:
   $$p_{\text{selective}} = 2 \min \left( \frac{1}{N} \sum_{i=1}^N \mathbb{I}(z_i \ge \phi_{\text{obs}}), \frac{1}{N} \sum_{i=1}^N \mathbb{I}(z_i \le \phi_{\text{obs}}) \right)$$

---

## 5. Implementation 2: K-Means Analytical Selective Inference

For K-Means, the selection event is defined by the final Voronoi cell assignments. A data point $X_i(z)$ belongs to cluster $k$ if and only if it is closer to the center of cluster $k$ than to any other center $j \neq k$:
$$\| X_i(z) - \mu_k(z) \|_2^2 \le \| X_i(z) - \mu_j(z) \|_2^2 \quad \forall j \neq k$$

Since both the perturbed points $X_i(z) = X_{\text{perp}, i} + z \nu_i$ and the centers $\mu(z)$ are linear functions of $z$, this boundary condition simplifies to a set of quadratic inequalities in $z$:
$$c_2 z^2 + c_1 z + c_0 \le 0$$
where:
$$c_2 = \|a_{\nu}\|_2^2 - \|b_{\nu}\|_2^2$$
$$c_1 = 2 ( \langle a_{\text{perp}}, a_{\nu} \rangle - \langle b_{\text{perp}}, b_{\nu} \rangle )$$
$$c_0 = \|a_{\text{perp}}\|_2^2 - \|b_{\text{perp}}\|_2^2$$
and:
- $a(z) = X_i(z) - \mu_k(z)$
- $b(z) = X_i(z) - \mu_j(z)$

In [kmeans_inference.py](file:///Users/halu/Code/tme_analysis/scripts/selective_inference/kmeans_inference.py#L164), solving these quadratic inequalities analytically across all points $i$ and alternative clusters $j$ yields the exact truncation interval:
$$\mathcal{S} = [V_{\text{min}}, V_{\text{max}}]$$

Because $\mathcal{S}$ is an interval, the conditional null distribution is a **Truncated Normal Distribution**:
$$z \mid \{C(X(z)) = \hat{y}_{\text{obs}}\} \;\sim\; \text{TruncatedNormal}(0, \sigma_z^2, V_{\text{min}}, V_{\text{max}})$$
The analytical selective $p$-value is then computed using the cumulative distribution function $F(z)$ of this truncated normal:
$$p_{\text{selective}} = 2 \min(F(z), 1 - F(z))$$
where:
$$F(z) = \frac{\Phi(z / \sigma_z) - \Phi(V_{\text{min}} / \sigma_z)}{\Phi(V_{\text{max}} / \sigma_z) - \Phi(V_{\text{min}} / \sigma_z)}$$
and $\Phi$ is the standard normal CDF.
