"""Core enumeration types and typeclasses for tme_response."""

from __future__ import annotations

from enum import Enum, unique


@unique
class PredictorCategory(str, Enum):
    """Categorization of ICB response predictors."""
    SIGNATURE = "signature"
    CELL_DECONVOLUTION = "cell_deconvolution"
    SYSTEMS_MODEL = "systems_model"
    GENOMIC = "genomic"
    COMPOSITE_SYNERGY = "composite_synergy"


@unique
class ValidationMetric(str, Enum):
    """Evaluation metrics for biomarker concordance and discrimination."""
    ROC_AUC = "roc_auc"
    PR_AUC = "pr_auc"
    BRIER_SCORE = "brier_score"
    ODDS_RATIO = "odds_ratio"
    CONCORDANCE_INDEX = "concordance_index"
    COX_HAZARD_RATIO = "cox_hazard_ratio"


@unique
class BiopsyTimepoint(str, Enum):
    """Standardized biopsy timepoint relative to treatment initiation."""
    PRE = "Pre"
    ON_TREATMENT = "On-Treatment"
    POST = "Post"
    UNKNOWN = "Unknown"


@unique
class ClinicalResponse(str, Enum):
    """Standardized RECIST/clinical benefit categorization."""
    RESPONDER = "Responder"
    NON_RESPONDER = "Non-Responder"
    UNKNOWN = "Unknown"
