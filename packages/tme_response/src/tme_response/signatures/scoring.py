"""Unified signature computation and consolidation engine."""

from __future__ import annotations

import polars as pl
from returns.result import Failure, Result, Success

from ..schemas import MultiOmicCohort
from .compass_baselines import (
    compute_ayers_ifng6_score,
    compute_cd8_duo_score,
    compute_cristescu_gep_score,
    compute_davoli_cis_score,
    compute_fehrenbacher_teff_score,
    compute_freeman_pgm_score,
    compute_genebio_target_score,
    compute_huang_nrs_score,
    compute_jiang_ctls_score,
    compute_jiang_tams_score,
    compute_jiang_texh_score,
    compute_kong_netbio_score,
    compute_messina_cks_score,
    compute_nurmik_cafs_score,
    compute_roh_is_score,
    compute_wu_mias_score,
)
from .cyt import compute_cyt_score
from .gep import compute_gep_score
from .impres import compute_impres_score
from .ipres import compute_ipres_score
from .single_gene import compute_single_gene_score


def compute_standard_signatures(cohort: MultiOmicCohort) -> Result[pl.DataFrame, str]:
    """Compute all standard transcriptomic signatures for a cohort and return an aligned Polars DataFrame.

    Output columns include standard signatures plus all COMPASS Table S2 baselines.
    """
    try:
        base_df = pl.DataFrame({"sample_id": list(cohort.sample_ids)})

        signature_funcs = [
            ("score_CYT", lambda: compute_cyt_score(cohort)),
            ("score_IMPRES", lambda: compute_impres_score(cohort)),
            ("score_Ayers_GEP", lambda: compute_gep_score(cohort)),
            ("score_CXCL9", lambda: compute_single_gene_score(cohort, "CXCL9")),
            ("score_CD8A", lambda: compute_single_gene_score(cohort, "CD8A")),
            ("score_PDCD1", lambda: compute_single_gene_score(cohort, "PDCD1")),
            ("score_CD274", lambda: compute_single_gene_score(cohort, "CD274")),
            ("score_CTLA4", lambda: compute_single_gene_score(cohort, "CTLA4")),
            ("score_IPRES", lambda: compute_ipres_score(cohort, invert_for_response=True)),
            ("score_GeneBio", lambda: compute_genebio_target_score(cohort)),
            ("score_CD8_Duo", lambda: compute_cd8_duo_score(cohort)),
            ("score_Davoli_CIS", lambda: compute_davoli_cis_score(cohort)),
            ("score_Fehrenbacher_Teff", lambda: compute_fehrenbacher_teff_score(cohort)),
            ("score_Freeman_PGM", lambda: compute_freeman_pgm_score(cohort)),
            ("score_Huang_NRS", lambda: compute_huang_nrs_score(cohort)),
            ("score_Ayers_IFNG_6", lambda: compute_ayers_ifng6_score(cohort)),
            ("score_Jiang_CTLs", lambda: compute_jiang_ctls_score(cohort)),
            ("score_Jiang_TAMs_Inverted", lambda: compute_jiang_tams_score(cohort, invert_for_response=True)),
            ("score_Jiang_Texh_Inverted", lambda: compute_jiang_texh_score(cohort, invert_for_response=True)),
            ("score_Messina_CKS", lambda: compute_messina_cks_score(cohort)),
            ("score_Nurmik_CAFs_Inverted", lambda: compute_nurmik_cafs_score(cohort, invert_for_response=True)),
            ("score_Roh_IS", lambda: compute_roh_is_score(cohort)),
            ("score_Wu_MIAS", lambda: compute_wu_mias_score(cohort)),
            ("score_Cristescu_GEP", lambda: compute_cristescu_gep_score(cohort)),
            ("score_Kong_NetBio_PD1", lambda: compute_kong_netbio_score(cohort)),
        ]

        for col_name, func in signature_funcs:
            match func():
                case Success(pred):
                    base_df = base_df.join(
                        pred.predictions.rename({"score": col_name}),
                        on="sample_id",
                        how="left",
                    )
                case Failure(_):
                    pass

        return Success(base_df)
    except Exception as exc:
        return Failure(f"Failed to consolidate signatures for {cohort.cohort_id}: {exc}")

