"""Malignant cell sampling and tumor compartment balance configuration."""

from __future__ import annotations

from enum import Enum
from typing import Mapping, Sequence
import anndata as ad  # type: ignore[import-untyped]
import numpy as np
import scipy.sparse as sp
from pydantic import BaseModel, ConfigDict, Field
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success

from ..deconvolution.malignant import calculate_cnv_proxy_scores, detect_malignant_cells
from ..deconvolution.models import DeconvolutionReferenceConfig
from ..logging import get_logger

logger = get_logger("sampling.malignant")


class MalignantStrategy(str, Enum):
    """Strategy for incorporating tumor/malignant cells into deconvolution reference."""

    POOLED_GENERIC = "pooled_generic"
    PATIENT_STRATIFIED = "patient_stratified"
    EXCLUDED_TME_ONLY = "excluded_tme_only"


class MalignantSamplingConfig(BaseModel):
    """Immutable specification for malignant vs non-malignant cell sampling."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    strategy: MalignantStrategy = MalignantStrategy.PATIENT_STRATIFIED
    malignant_fraction: float = 0.20
    max_cells_per_patient: int = 50
    contig_variance_qc: bool = True
    tumor_marker_threshold: float = 0.50
    patient_id_keys: tuple[str, ...] = ("patient_id", "patient", "donor_id", "donor", "sample_id", "orig.ident")


def sample_malignant_and_tme_cells(
    adata: ad.AnnData,
    config: MalignantSamplingConfig,
    target_total_cells: int,
    ref_config: DeconvolutionReferenceConfig = DeconvolutionReferenceConfig(),
    seed: int = 42,
) -> Result[ad.AnnData, str]:
    """Sample cells balancing malignant and microenvironmental (immune/stromal) compartments.

    Args:
        adata: Input AnnData dataset.
        config: MalignantSamplingConfig defining strategy, malignant fraction, and patient cap.
        target_total_cells: Total cell budget to draw from this dataset.
        ref_config: DeconvolutionReferenceConfig for malignant detection heuristics.
        seed: Random seed for reproducibility.

    Returns:
        Success(subsampled_adata) or Failure(error message).
    """
    if adata.n_obs == 0:
        return Failure("Cannot sample from an empty AnnData object.")

    rng = np.random.default_rng(seed)

    # 1. Identify malignant cells
    is_mal = detect_malignant_cells(adata, ref_config)
    n_mal_available = int(np.sum(is_mal))
    n_tme_available = int(adata.n_obs - n_mal_available)

    logger.info(
        "Malignant detection: %d malignant cells, %d non-malignant TME cells available",
        n_mal_available,
        n_tme_available,
    )

    all_indices = np.arange(adata.n_obs)
    mal_indices = all_indices[is_mal]
    tme_indices = all_indices[~is_mal]

    # Handle CNV contig variance QC filtering if enabled
    if config.contig_variance_qc and n_mal_available > 10:
        cnv_scores = calculate_cnv_proxy_scores(adata)
        if cnv_scores is not None:
            # Retain malignant cells in top 80% CNV score to purge misclassified normal cells
            mal_cnv = cnv_scores[is_mal]
            q20 = float(np.percentile(mal_cnv, 20))
            valid_mal_mask = mal_cnv >= q20
            mal_indices = mal_indices[valid_mal_mask]
            n_mal_available = len(mal_indices)

    # 2. Dispatch based on strategy
    match config.strategy:
        case MalignantStrategy.EXCLUDED_TME_ONLY:
            # Zero malignant cells; sample entirely from TME
            if n_tme_available == 0:
                return Failure("Strategy is EXCLUDED_TME_ONLY but no non-malignant cells detected.")
            n_to_sample = min(target_total_cells, n_tme_available)
            chosen_tme = rng.choice(tme_indices, size=n_to_sample, replace=False)
            chosen_all = sorted(chosen_tme)

        case MalignantStrategy.POOLED_GENERIC:
            # Target specified malignant fraction
            target_mal = int(round(target_total_cells * config.malignant_fraction))
            n_mal = max(0, min(target_mal, n_mal_available))
            n_tme = min(target_total_cells - n_mal, n_tme_available)

            chosen_mal = rng.choice(mal_indices, size=n_mal, replace=False) if n_mal > 0 else np.array([], dtype=int)
            chosen_tme = rng.choice(tme_indices, size=n_tme, replace=False) if n_tme > 0 else np.array([], dtype=int)
            chosen_all = sorted(np.concatenate([chosen_mal, chosen_tme]).astype(int))

        case MalignantStrategy.PATIENT_STRATIFIED:
            # Identify patient / donor metadata column
            patient_col = next((c for c in config.patient_id_keys if c in adata.obs.columns), None)

            target_mal = int(round(target_total_cells * config.malignant_fraction))

            if patient_col is not None and n_mal_available > 0 and target_mal > 0:
                patient_vals = np.asarray(adata.obs[patient_col])[mal_indices]
                unique_patients = np.unique(patient_vals)

                selected_mal_list: list[int] = []
                # Allocate cells per patient
                per_patient_quota = max(1, min(config.max_cells_per_patient, target_mal // max(1, len(unique_patients))))

                for pat in unique_patients:
                    pat_mask = patient_vals == pat
                    pat_global_indices = mal_indices[pat_mask]
                    n_draw = min(len(pat_global_indices), per_patient_quota)
                    if n_draw > 0:
                        sampled_pat = rng.choice(pat_global_indices, size=n_draw, replace=False)
                        selected_mal_list.extend(sampled_pat)

                # If still below target_mal, draw remaining uniformly
                if len(selected_mal_list) < target_mal:
                    remaining_pool = list(set(mal_indices) - set(selected_mal_list))
                    extra_needed = min(target_mal - len(selected_mal_list), len(remaining_pool))
                    if extra_needed > 0:
                        extra_drawn = rng.choice(remaining_pool, size=extra_needed, replace=False)
                        selected_mal_list.extend(extra_drawn)

                chosen_mal = np.array(selected_mal_list[:target_mal], dtype=int)
            else:
                n_mal = max(0, min(target_mal, n_mal_available))
                chosen_mal = rng.choice(mal_indices, size=n_mal, replace=False) if n_mal > 0 else np.array([], dtype=int)

            n_tme = min(target_total_cells - len(chosen_mal), n_tme_available)
            chosen_tme = rng.choice(tme_indices, size=n_tme, replace=False) if n_tme > 0 else np.array([], dtype=int)
            chosen_all = sorted(np.concatenate([chosen_mal, chosen_tme]).astype(int))

    if not chosen_all:
        return Failure("No cells could be selected under the specified malignant sampling configuration.")

    sub_adata = adata[chosen_all].copy()
    sub_adata.obs["is_malignant"] = is_mal[chosen_all]
    return Success(sub_adata)
