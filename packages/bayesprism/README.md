# BayesPrism (Python)

A pure functional Python and PyTorch implementation of **BayesPrism** for Bayesian cell type and gene expression deconvolution using single-cell RNA-seq references.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

---

## Overview

BayesPrism infers the cellular composition and cell-type-specific gene expression profiles of heterogeneous bulk RNA-seq data by leveraging single-cell transcriptomics as a prior. 

This package is a high-performance Python port built upon:
- **PyTorch** for accelerated matrix computations and GPU tensor operations.
- **Polars** for fast, declarative data processing.
- **Pure Functional Design** ensuring immutable state transitions, railway-oriented error handling via `returns`, and strict type safety.

---

## Installation

```bash
pip install bayesprism
```

Or install directly from GitHub:

```bash
pip install git+https://github.com/LuJoHae/bayesprism-python.git
```

---

## Quick Start

```python
import numpy as np
from bayesprism import new_prism, run_prism, get_fraction, get_exp

# 1. Prepare raw count matrices (samples x genes, cell_types x genes)
# Bulk mixture: X (samples x genes)
# scRNA-seq reference: phi_ref (cell_types x genes)
gene_names = ("CD3D", "CD8A", "CD4", "MS4A1", "CD68", ...)

# 2. Initialize the Prism model container
prism_result = new_prism(
    X_bulk=X_bulk_counts,
    phi_ref=phi_sc_reference,
    gene_names=gene_names,
)

# Unwrap monadic Result safely
prism = prism_result.unwrap()

# 3. Run Bayesian MCMC deconvolution
fitted_prism = run_prism(prism)

# 4. Extract inferred cell type fractions (theta)
fractions = get_fraction(fitted_prism)

# 5. Extract cell-type-specific expression profiles (Z)
expression_profiles = get_exp(fitted_prism)
```

---

## Key Features

- **Bayesian Joint Estimation**: Jointly infers cell-type proportions ($\theta$) and cell-type-specific gene expression ($Z$).
- **Linear Count Space**: Operates directly on unnormalized raw sequencing counts to preserve count distributions without log-transformation artifacts.
- **Batched GPU Acceleration**: PyTorch tensors enable fast vectorized updates across large patient cohorts.
- **Modular Pipeline**: Full functional pipeline with separate Gibbs sampling, Nelder-Mead reference updating, and NMF embedding learning.

---

## Citation & Acknowledgments

If you use BayesPrism in your research, please cite the original BayesPrism methodology:

> Chu, T., Wang, Z., Pe'er, D., & Danko, C. G. (2022). Cell type and gene expression deconvolution with BayesPrism. *Nature Cancer*, 3(4), 505-517. [doi:10.1038/s43018-022-00356-3](https://doi.org/10.1038/s43018-022-00356-3)

---

## License

This project is licensed under the [MIT License](LICENSE).
