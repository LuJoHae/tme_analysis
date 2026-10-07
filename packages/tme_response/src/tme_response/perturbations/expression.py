"""Pure functional transcriptomic perturbation operators on MultiOmicCohort models."""

from __future__ import annotations

from typing import Sequence
import numpy as np
import polars as pl
from tme_datasets.genesets import AYERS_T_CELL_INFLAMED_GEP

from ..schemas import MultiOmicCohort
from .models import DilutionConfig, DropoutConfig, JitterConfig

DEFAULT_IMMUNE_EFFECTOR_GENES = tuple(set(
    tuple(AYERS_T_CELL_INFLAMED_GEP.genes)
    + ("GZMA", "PRF1", "CXCL9", "CD8A", "PDCD1", "IFNG", "CTLA4", "LAG3", "HAVCR2")
))


def apply_expression_jitter(
    cohort: MultiOmicCohort,
    config: JitterConfig,
) -> MultiOmicCohort:
    """Apply multiplicative log-normal expression jitter to non-zero TPM values.

    TPM' = TPM * 2^(N(0, sigma^2)). Simulates technical assay variance and quantification noise.
    """
    if config.sigma <= 0.0:
        return cohort

    rng = np.random.default_rng(config.seed)
    gene_cols = [c for c in cohort.expression_tpm.columns if c != "sample_id"]
    if not gene_cols:
        return cohort

    mat = cohort.expression_tpm.select(gene_cols).to_numpy()
    noise = 2.0 ** rng.normal(0.0, config.sigma, size=mat.shape)
    perturbed_mat = np.maximum(mat * noise, 0.0)

    col_dict: dict[str, list[object]] = {
        "sample_id": cohort.expression_tpm["sample_id"].to_list()
    }
    for idx, g in enumerate(gene_cols):
        col_dict[g] = [float(v) for v in perturbed_mat[:, idx]]

    perturbed_df = pl.DataFrame(col_dict)
    return cohort.model_copy(update={"expression_tpm": perturbed_df})


def apply_gene_dropout(
    cohort: MultiOmicCohort,
    config: DropoutConfig,
) -> MultiOmicCohort:
    """Simulate technical gene dropout by setting values to 0 with probability dropout_rate.

    Simulates shallow coverage, probe hybridization dropouts, or targeted panel gene omission.
    """
    if config.dropout_rate <= 0.0:
        return cohort

    rng = np.random.default_rng(config.seed)
    gene_cols = [c for c in cohort.expression_tpm.columns if c != "sample_id"]
    if not gene_cols:
        return cohort

    mat = cohort.expression_tpm.select(gene_cols).to_numpy()
    rate = min(1.0, max(0.0, config.dropout_rate))
    dropout_mask = rng.uniform(0.0, 1.0, size=mat.shape) < rate
    perturbed_mat = np.where(dropout_mask, 0.0, mat)

    col_dict: dict[str, list[object]] = {
        "sample_id": cohort.expression_tpm["sample_id"].to_list()
    }
    for idx, g in enumerate(gene_cols):
        col_dict[g] = [float(v) for v in perturbed_mat[:, idx]]

    perturbed_df = pl.DataFrame(col_dict)
    return cohort.model_copy(update={"expression_tpm": perturbed_df})


def apply_immune_dilution(
    cohort: MultiOmicCohort,
    config: DilutionConfig,
    immune_genes: Sequence[str] | None = None,
) -> MultiOmicCohort:
    """Simulate stromal infiltration and low tumor purity by downscaling the immune compartment.

    TPM_immune' = dilution_factor * TPM_immune.
    dilution_factor = 1.0 indicates pure baseline; 0.2 indicates severe (80%) stromal dilution.
    """
    if config.dilution_factor >= 1.0:
        return cohort

    alpha = float(max(0.0, min(1.0, config.dilution_factor)))
    target_immune = set(immune_genes or DEFAULT_IMMUNE_EFFECTOR_GENES)

    gene_cols = [c for c in cohort.expression_tpm.columns if c != "sample_id"]
    if not gene_cols:
        return cohort

    mat = cohort.expression_tpm.select(gene_cols).to_numpy()
    scale_vector = np.array(
        [alpha if g in target_immune else 1.0 for g in gene_cols],
        dtype=np.float32,
    )
    perturbed_mat = mat * scale_vector.reshape(1, -1)

    col_dict: dict[str, list[object]] = {
        "sample_id": cohort.expression_tpm["sample_id"].to_list()
    }
    for idx, g in enumerate(gene_cols):
        col_dict[g] = [float(v) for v in perturbed_mat[:, idx]]

    perturbed_df = pl.DataFrame(col_dict)
    return cohort.model_copy(update={"expression_tpm": perturbed_df})
