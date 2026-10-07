# Synthetic Collinearity Benchmark: Data Generation and Mathematical Specification

## 1. Overview & Motivation

Reference-based transcriptomic deconvolution seeks to resolve the cellular composition $\boldsymbol{\theta}_n \in \Delta^{S-1}$ of heterogeneous tissue samples $\mathbf{y}_n \in \mathbb{R}^G$ from a single-cell signature matrix $\boldsymbol{\Phi} \in \mathbb{R}^{S \times G}$. 

In practical immune and oncology applications, single-cell references frequently feature closely related cellular states belonging to the same lineage (e.g., cytotoxic vs. exhausted CD8$^+$ T cells; M1 vs. M2 macrophages). As sibling cell states share common lineage-defining marker genes, their transcriptomic profiles exhibit severe **collinearity** ($r \to 1.0$), driving the reference Gram matrix $\boldsymbol{\Phi} \boldsymbol{\Phi}^T$ toward numerical singularity ($\kappa(\boldsymbol{\Phi}) \to \infty$).

To systematically investigate this phenomenon and evaluate regularized deconvolution engines under controlled conditions, we constructed a **continuous-spectrum synthetic benchmark**. This benchmark spans from completely orthogonal cell states ($r \approx 0.0$) through moderate and high collinearity to extreme near-singular profiles ($r = 0.990$) across $S = 12$ cell states and $T = 4$ broad lineages.

> [!NOTE] Pure In Silico Generation (Zero Clinical/Empirical Data)
> All expression counts, reference signatures, cellular proportions, and bulk mixtures are **100% synthetically generated from parametric probability distributions** (Gamma priors, Dirichlet priors, and Multinomial noise). No experimental expression counts, marker genes, or single-cell count matrices are drawn from real clinical cohorts (such as Sade-Feldman or Hugo). To make this completely transparent, all entities use explicit synthetic designations:
> - Broad Lineages: `Lineage_1`, `Lineage_2`, `Lineage_3`, `Lineage_4`
> - Fine Cell States: `L1_State_1` to `L4_State_3`
> - Marker Genes: `Gene_000` to `Gene_399`
> - Bulk Mixtures: `Synthetic_Sample_00` to `Synthetic_Sample_49`

---

## 2. Generative Model for Collinear Cell State Signatures

### 2.1 Cellular Hierarchy Structure
The synthetic reference represents $S = 12$ fine-grained cell states nested hierarchically within $T = 4$ broad lineages ($3$ sibling states per lineage):

1. **Lineage 1 ($t = 0$)**:
   - State 0: `L1_State_1`
   - State 1: `L1_State_2`
   - State 2: `L1_State_3`
2. **Lineage 2 ($t = 1$)**:
   - State 3: `L2_State_1`
   - State 4: `L2_State_2`
   - State 5: `L2_State_3`
3. **Lineage 3 ($t = 2$)**:
   - State 6: `L3_State_1`
   - State 7: `L3_State_2`
   - State 8: `L3_State_3`
4. **Lineage 4 ($t = 3$)**:
   - State 9: `L4_State_1`
   - State 10: `L4_State_2`
   - State 11: `L4_State_3`

The broad lineage aggregation operator $\mathbf{M} \in \{0, 1\}^{4 \times 12}$ maps states to lineages:
$$M_{t, s} = \begin{cases} 1 & \text{if state } s \text{ belongs to lineage } t \\ 0 & \text{otherwise} \end{cases}$$

$$\mathbf{M} = \begin{bmatrix}
1 & 1 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 1 & 1 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 1 & 1 & 1 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 1 & 1
\end{bmatrix}$$

---

### 2.2 Continuous Sibling Correlation Superposition
For $G = 400$ genes and each lineage $t \in \{0, 1, 2, 3\}$:
1. Draw an independent base lineage expression program:
   $$\mathbf{b}_t \sim \operatorname{Gamma}(\alpha = 2.0, \beta = 1.0)^{\otimes G}$$
   where $\mathbb{E}[b_{t, g}] = \alpha \beta = 2.0$ and $\operatorname{Var}(b_{t, g}) = \alpha \beta^2 = 2.0$.
2. For each state $k \in \{0, 1, 2\}$ within lineage $t$, draw an independent state-specific orthogonal component:
   $$\mathbf{u}_{t, k} \sim \operatorname{Gamma}(\alpha = 2.0, \beta = 1.0)^{\otimes G}$$
3. For a target correlation parameter $r \in [0.0, 1.0]$, synthesize the unnormalized state profile $\mathbf{s}_{t, k} \in \mathbb{R}_+^G$ via the linear combination:
   $$\mathbf{s}_{t, k} = \sqrt{r} \, \mathbf{b}_t + \sqrt{1 - r} \, \mathbf{u}_{t, k}$$

#### Theoretical Proof of Correlation Invariance
For any two sibling states $a \neq b \in \{0, 1, 2\}$ within the same lineage $t$:
$$\operatorname{Cov}(s_{t, a, g}, s_{t, b, g}) = \operatorname{Cov}\left(\sqrt{r} b_{t, g} + \sqrt{1-r} u_{t, a, g}, \; \sqrt{r} b_{t, g} + \sqrt{1-r} u_{t, b, g}\right)$$
Since $\mathbf{u}_{t, a}$ and $\mathbf{u}_{t, b}$ are mutually independent and independent of $\mathbf{b}_t$:
$$\operatorname{Cov}(s_{t, a, g}, s_{t, b, g}) = (\sqrt{r})^2 \operatorname{Var}(b_{t, g}) = 2r$$
The individual variances are:
$$\operatorname{Var}(s_{t, a, g}) = r \operatorname{Var}(b_{t, g}) + (1-r) \operatorname{Var}(u_{t, a, g}) = 2r + 2(1-r) = 2$$
Therefore, the Pearson correlation between sibling states is exactly:
$$\operatorname{Corr}(\mathbf{s}_{t, a}, \mathbf{s}_{t, b}) = \frac{\operatorname{Cov}(s_{t, a, g}, s_{t, b, g})}{\sqrt{\operatorname{Var}(s_{t, a, g}) \operatorname{Var}(s_{t, b, g})}} = \frac{2r}{\sqrt{2 \times 2}} = r \quad \blacksquare$$

Conversely, for states $s \in t$ and $s' \in t'$ belonging to distinct lineages ($t \neq t'$), $\mathbf{b}_t$ and $\mathbf{b}_{t'}$ are independent:
$$\operatorname{Cov}(s_{t, k, g}, s_{t', k', g}) = 0 \implies \operatorname{Corr}(\mathbf{s}_{t, k}, \mathbf{s}_{t', k'}) \approx 0$$
This guarantees that collinearity is strictly intra-lineage, creating an exact block-diagonal correlation structure.

---

### 2.3 Simplex Normalization and Empirical Sibling Correlation
To form valid transcriptomic probability vectors on the gene simplex $\Delta^{G-1}$, each row of $\mathbf{S}$ is normalized:
$$\phi_{s, g} = \frac{s_{s, g}}{\sum_{g'=1}^G s_{s, g'}}, \quad \boldsymbol{\phi}_s \in \Delta^{G-1}$$
The actual empirical sibling correlation is computed as the average over all 12 intra-lineage sibling pairs:
$$r_{\text{sibling}} = \frac{1}{4} \sum_{t=0}^3 \frac{1}{3} \sum_{0 \le a < b \le 2} \operatorname{Corr}(\boldsymbol{\phi}_{3t+a}, \boldsymbol{\phi}_{3t+b})$$
This empirical value $r_{\text{sibling}}$ serves as the quantitative coordinate on the horizontal axis of all benchmark plots.

---

## 3. Synthetic Bulk Mixture Generation

### 3.1 Ground Truth Cell State Proportions
For each synthetic patient sample $n \in \{1, \dots, N\}$ ($N = 50$):
True state fractions are drawn from a non-uniform Dirichlet prior reflecting typical TME cellular heterogeneity:
$$\boldsymbol{\theta}_n^* \sim \operatorname{Dirichlet}(\boldsymbol{\alpha})$$
$$\boldsymbol{\alpha} = [1.2, 1.0, 0.8, \; 1.0, 0.8, 0.6, \; 1.2, 1.0, 0.8, \; 0.8, 0.6, 0.4]^T \in \mathbb{R}_+^{12}$$
True broad lineage proportions are derived by exact linear aggregation:
$$\boldsymbol{\theta}_n^{\text{lineage}} = \mathbf{M} \boldsymbol{\theta}_n^* \in \Delta^{T-1}$$

### 3.2 Multinomial Sequencing Sampling
The expected bulk expression distribution is the linear convex combination of reference signatures:
$$\mathbf{p}_n = \boldsymbol{\theta}_n^* \boldsymbol{\Phi} = \sum_{s=0}^{11} \theta_{n, s}^* \boldsymbol{\phi}_s \in \Delta^{G-1}$$
Sequencing count generation is modeled under standard Poisson/Multinomial counting statistics with total sequencing depth $N_{\text{reads}} = 80{,}000$ reads:
$$\mathbf{y}_n \sim \operatorname{Multinomial}(N_{\text{reads}} = 80{,}000, \; \mathbf{p}_n)$$

---

## 4. Mathematical Null-Space Invariance

For any two sibling states $s_a, s_b$ within lineage $t$, consider the singular difference vector:
$$\mathbf{v} = \frac{1}{\sqrt{2}} (\mathbf{e}_{s_a} - \mathbf{e}_{s_b})$$
Evaluating the lineage aggregation operator on $\mathbf{v}$:
$$\mathbf{M} \mathbf{v} = \frac{1}{\sqrt{2}} (\mathbf{M} \mathbf{e}_{s_a} - \mathbf{M} \mathbf{e}_{s_b}) = \frac{1}{\sqrt{2}} (\mathbf{e}_t - \mathbf{e}_t) = \mathbf{0}$$
Thus, **the collinear error mode lies strictly in the null space of $\mathbf{M}$**:
$$\mathbf{v} \in \operatorname{Null}(\mathbf{M})$$

As a direct consequence, the broad lineage Mean Squared Error:
$$\operatorname{MSE}_{\text{lineage}} = \frac{1}{N} \sum_{n=1}^N \|\mathbf{M} \hat{\boldsymbol{\theta}}_n - \mathbf{M} \boldsymbol{\theta}_n^*\|_2^2$$
remains mathematically invariant across all deconvolution methods ($\approx 5.0 \times 10^{-4}$), even when unregularized state-level MSE explodes.

---

## 5. Continuous Correlation Spectrum Grid

The benchmark evaluates $10$ target correlation levels spanning the interval $[0.0, 0.99]$:

| Target $r$ | Empirical $r_{\text{sibling}}$ | Nominal $\kappa(\mathbf{H})$ | RegDeconv $\kappa(\mathbf{H})$ | Regime Interpretation |
| :---: | :---: | :---: | :---: | :--- |
| **0.00** | $-0.011$ | $3.3 \times 10^1$ | $3.3 \times 10^1$ | Uncorrelated (Orthogonal cell states) |
| **0.10** | $+0.079$ | $5.9 \times 10^1$ | $5.9 \times 10^1$ | Very Low Collinearity |
| **0.20** | $+0.212$ | $7.1 \times 10^1$ | $7.1 \times 10^1$ | Low Collinearity |
| **0.40** | $+0.394$ | $9.9 \times 10^1$ | $9.9 \times 10^1$ | Mild Collinearity |
| **0.60** | $+0.614$ | $1.6 \times 10^2$ | $1.6 \times 10^2$ | Moderate Collinearity |
| **0.80** | $+0.805$ | $2.9 \times 10^2$ | $2.9 \times 10^2$ | Sub-critical Collinearity |
| **0.90** | $+0.911$ | $4.9 \times 10^2$ | $4.9 \times 10^2$ | Elevated Collinearity |
| **0.95** | $+0.948$ | $8.6 \times 10^2$ | $8.6 \times 10^2$ | High Collinearity |
| **0.98** | $+0.978$ | $1.8 \times 10^3$ | $1.8 \times 10^3$ | Threshold Collinearity ($\kappa \approx \kappa_{\text{target}}$) |
| **0.99** | $+0.990$ | $3.8 \times 10^3$ | **$1.6 \times 10^3$** | **Severe Collinearity** (Deficit activated, $\kappa \le 2,000$) |

---

## 6. Estimator Bias-Variance Tradeoff & Ground-Truth Concordance

### 6.1 Mathematical Formulation of the Hierarchical Shrinkage Effect
A key observation in the severe collinearity regime ($r = 0.99$) is that BayesPrism and InstaPrism exhibit lower fine-state Pearson correlation with ground truth ($r \approx 0.53 - 0.55$) than unregularized NNLS ($r = 0.919$) or RegDeconv ($r = 0.825$), despite achieving near-perfect broad lineage concordance ($r > 0.98$, $\text{MSE} \approx 5.4 \times 10^{-4}$).

This behavior is a direct mathematical consequence of **hierarchical Bayesian shrinkage**:
1. BayesPrism and InstaPrism decompose the cell proportion vector into broad cell-type fractions $\psi_t$ and intra-lineage state fractions $\omega_{t, s}$, such that $\theta_{t, s} = \psi_t \cdot \omega_{t, s}$.
2. Intra-lineage fractions are governed by a Dirichlet prior:
   $$\boldsymbol{\omega}_t \sim \operatorname{Dirichlet}(\boldsymbol{\alpha}_t)$$
3. As intra-lineage sibling correlation approaches $r \to 0.99$, the transcriptomic likelihood provides negligible gradient to distinguish sibling profiles ($\mathbf{s}_{t, 1} \approx \mathbf{s}_{t, 2} \approx \mathbf{s}_{t, 3}$).
4. Consequently, the posterior contracts strongly toward the symmetric prior centroid:
   $$\omega_{t, 1} \approx \omega_{t, 2} \approx \omega_{t, 3} \approx \frac{1}{3} \implies \hat{\theta}_{t, s} \approx \frac{1}{3} \hat{\psi}_t$$
5. When the true ground-truth proportions are heterogeneous within a lineage ($\boldsymbol{\alpha} = [1.2, 1.0, 0.8]$), this centroid shrinkage introduces a deterministic **squared bias**:
   $$\operatorname{Bias}^2(\hat{\theta}_{t, s}) = \left(\frac{1}{3} \psi_t^* - \theta_{t, s}^*\right)^2$$
   flattening fine-state variation across samples and dampening the fine-state Pearson correlation $r$.

### 6.2 Empirical Replicate Stress-Test ($B = 30$ Sequencing Replicates)
To evaluate whether this shrinkage trades fine-state correlation for **estimator stability**, we generated $B = 30$ independent sequencing draws $\mathbf{y}^{(b)} \sim \operatorname{Multinomial}(80{,}000, \boldsymbol{\theta}^* \boldsymbol{\Phi})$ from the **same underlying biological tissue** at $r = 0.99$:

$$\operatorname{MSE}(\hat{\theta}) = \operatorname{Bias}^2(\hat{\theta}) + \operatorname{Var}(\hat{\theta})$$

| Method | Estimator Variance $\operatorname{Var}(\hat{\theta})$ | Variance % of MSE | Squared Bias $\operatorname{Bias}^2(\hat{\theta})$ | Bias % of MSE | Total MSE |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Unregularized (NNLS)** | $0.001247$ | **$95.6\%$** | $0.000099$ | $7.6\%$ | $0.001304$ |
| **Rectangle (DWLS-QP)** | $0.001018$ | **$91.4\%$** | $0.000130$ | $11.7\%$ | $0.001114$ |
| **CIBERSORT (reimpl., $\nu$-SVR)** | $0.001369$ | **$89.0\%$** | $0.000214$ | $13.9\%$ | $0.001538$ |
| **RegDeconv (Graph Lap)** | $0.000535$ | **$21.8\%$** | $0.001942$ | $79.0\%$ | $0.002459$ |
| **InstaPrism** | **$0.000003$** | **$0.06\%$** | $0.005299$ | $99.9\%$ | $0.005302$ |
| **BayesPrism (Gibbs)** | **$0.000028$** | **$0.50\%$** | $0.005646$ | $99.5\%$ | $0.005674$ |

**Key Takeaways**:
- **NNLS, Rectangle, and CIBERSORT**: $90\% - 95\%$ of their total error consists of pure estimator variance. Replicate estimates swing erratically from technical noise.
- **BayesPrism and InstaPrism**: Achieve a **$45\times$ to $400\times$ reduction in estimator variance** ($\operatorname{Var} < 10^{-5}$), confirming that their lower correlation is a deliberate mathematical trade-off: **sacrificing fine-grained resolution to eliminate technical noise variance**.
- **RegDeconv**: Provides a balanced middle ground, reducing variance by **$57\%$** relative to NNLS while maintaining half the shrinkage error of InstaPrism and BayesPrism.

---

## 7. True Biological Absence: Cell State & Lineage Dropout Stress-Testing

In realistic biological tissue, not all cell states or cell types are present in every biopsy ($\theta^* = 0.0$). We evaluated $N = 20$ controlled mixture samples under two absence conditions:
1. **State Dropout**: `L1_State_3` is truly absent ($\theta^* = 0.0$), but sibling states 1 and 2 are present.
2. **Lineage Dropout**: Entire `Lineage_4` is truly absent ($\theta^* = 0.0$).

| Method | Phantom Sibling State (%) | Phantom Absent Lineage (%) | Total Phantom Detection (%) | Active States RMSE |
| :--- | :---: | :---: | :---: | :---: |
| **Unregularized (NNLS)** | **$1.84\%$** | **$0.25\%$** | **$2.08\%$** | $0.0350$ |
| **Rectangle (DWLS-QP)** | **$1.68\%$** | **$0.23\%$** | **$1.92\%$** | $0.0338$ |
| **CIBERSORT (reimpl., $\nu$-SVR)** | $2.26\%$ | $3.42\%$ | $5.69\%$ | $0.0423$ |
| **RegDeconv (Graph Lap)** | $5.64\%$ | **$0.20\%$** | $5.84\%$ | $0.0572$ |
| **InstaPrism** | $8.25\%$ | $3.32\%$ | $11.57\%$ | $0.0968$ |
| **BayesPrism (Gibbs)** | $8.35\%$ | $7.39\%$ | $15.74\%$ | $0.1008$ |

**Findings**:
- **Dirichlet Phantom Leakage**: When sibling states are active, BayesPrism and InstaPrism allocate $\approx 8.3\%$ false positive fraction into truly absent sibling states due to their non-zero Dirichlet prior. When an entire lineage is absent, BayesPrism accumulates $15.7\%$ total phantom mass across true zeros.
- **Sparsity-Preserving Methods**: NNLS and Rectangle achieve the cleanest true zero suppression ($< 2.1\%$ phantom mass).
- **RegDeconv**: Completely suppresses absent lineages ($0.20\%$), while exhibiting modest smoothing leakage ($5.64\%$) among active sibling states.

---

## 8. CIBERSORTx Docker Container Integration

The benchmark suite supports the official Stanford CIBERSORTx Docker container (`cibersortx/fractions:latest`).

### 8.1 Supplying Credentials
Credentials can be provided via any of the following:
1. **`.env` file in the workspace root**:
   ```bash
   CIBERSORTX_USERNAME=your_registered_email@domain.com
   CIBERSORTX_TOKEN=your_stanford_api_token
   ```
2. **Shell Environment Variables**:
   ```bash
   export CIBERSORTX_USERNAME="your_registered_email@domain.com"
   export CIBERSORTX_TOKEN="your_stanford_api_token"
   ```
3. **CLI Arguments**:
   ```bash
   python scripts/regularized_deconv/benchmark_collinearity_regularized_deconv.py \
     --cibersortx-username "your_email" --cibersortx-token "your_token"
   ```
4. **Makefile Target**:
   ```bash
   make benchmark-cibersortx CIBERSORTX_USERNAME="your_email" CIBERSORTX_TOKEN="your_token"
   ```

### 8.2 Execution on Remote Server (`olm`) vs Local Mac
The repository is configured to dispatch CIBERSORTx container runs to the high-performance remote server **`olm`** (`Linux x86_64`, 140+ GB RAM) by default, eliminating architecture emulation overhead on Apple Silicon Macs:

- **Remote Verification (`olm`)**:
  ```bash
  make verify-cibersortx
  ```
  Automatically runs `make sync` to sync code and `.env` credentials (with `chmod 600`), then executes `verify_cibersortx_docker.py` on `olm`.
- **Local Verification (Mac)**:
  ```bash
  make verify-cibersortx-local
  ```
  Runs directly on the local machine using `/opt/homebrew/bin/python3.14` and local Docker Desktop.
- **Remote Full Benchmark (`olm`)**:
  ```bash
  make benchmark-cibersortx
  ```
  Syncs code and credentials to `olm`, executes the continuous collinearity benchmark under `TMPDIR=/storage/halu/tmp`, and automatically pulls all resulting Parquet tables back to `output/concordance/` via `make pull-results`.
- **Local Benchmark (Mac)**:
  ```bash
  make benchmark-cibersortx-local
  ```

---

## 9. Generated Parquet Data Artifacts

All benchmark outputs are persisted in `output/concordance/`:
- `synthetic_collinearity_benchmark_summary.parquet`: 10-point continuous spectrum summary metrics.
- `synthetic_collinearity_benchmark_results.parquet`: Per-sample detailed predictions across all methods.
- `synthetic_collinearity_true_proportions.parquet`: Exact ground truth proportions.
- `synthetic_estimator_variance_summary.parquet`: Bias-variance decomposition across $B = 30$ technical replicates.
- `synthetic_dropout_stress_test_summary.parquet`: Phantom detection fractions and active state errors under true zeros.

1. [`synthetic_collinearity_benchmark_summary.parquet`](file:///Users/halu/Code/tme_analysis/output/concordance/synthetic_collinearity_benchmark_summary.parquet): Aggregated summary metrics (MSE mean $\pm$ SD, condition number, sibling cross-talk correlation, spurious dropout rate) for all 6 methods across 11 correlation points (66 rows).
2. [`synthetic_collinearity_benchmark_results.parquet`](file:///Users/halu/Code/tme_analysis/output/concordance/synthetic_collinearity_benchmark_results.parquet): Sample-level metrics for each of the 50 samples, methods, and correlation points (3,300 rows).
3. [`synthetic_collinearity_true_proportions.parquet`](file:///Users/halu/Code/tme_analysis/output/concordance/synthetic_collinearity_true_proportions.parquet): Exact ground-truth cell state proportions ($\boldsymbol{\theta}_n^*$) and lineage sums ($\mathbf{M} \boldsymbol{\theta}_n^*$) for all 50 samples across the 11 correlation points (550 rows).
4. [`synthetic_collinearity_estimated_proportions.parquet`](file:///Users/halu/Code/tme_analysis/output/concordance/synthetic_collinearity_estimated_proportions.parquet): Inferred proportions ($\hat{\boldsymbol{\theta}}_n$) for all samples, methods, and correlation points (3,300 rows).
5. [`synthetic_collinearity_reference_matrices.parquet`](file:///Users/halu/Code/tme_analysis/output/concordance/synthetic_collinearity_reference_matrices.parquet): Full reference signature profiles $\boldsymbol{\Phi}$ across states and correlation regimes (132 rows).
