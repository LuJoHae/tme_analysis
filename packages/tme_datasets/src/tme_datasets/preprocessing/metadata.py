"""Metadata harmonization for clinical response and biological annotations."""

from __future__ import annotations

import numpy as np
import pandas as pd
import polars as pl
from returns.result import Failure, Result, Success

RESPONDER_TOKENS = frozenset({
    "CR", "PR", "RESPONDER", "RESPONSE", "YES", "TRUE", "1", "1.0", "DCB", "R",
    "COMPLETE RESPONSE", "PARTIAL RESPONSE",
})
NON_RESPONDER_TOKENS = frozenset({
    "PD", "SD", "NON-RESPONDER", "NON_RESPONDER", "NR", "NO", "FALSE", "0", "0.0", "NDB",
    "PROGRESSIVE DISEASE", "STABLE DISEASE", "PROGRESSION",
})


def binarize_response(val: object) -> float:
    """Standardize heterogeneous clinical annotations into binary 1.0 (Responder) vs 0.0 (Non-Responder)."""
    if val is None:
        return np.nan
    if isinstance(val, (int, float, np.integer, np.floating)):
        if np.isnan(val):
            return np.nan
        if val == 1.0:
            return 1.0
        if val == 0.0:
            return 0.0
    clean = str(val).strip().upper()
    if clean in RESPONDER_TOKENS:
        return 1.0
    if clean in NON_RESPONDER_TOKENS:
        return 0.0
    return np.nan


def standardize_recist(val: object) -> str:
    """Map string annotation into standard RECIST 1.1 category."""
    if val is None:
        return "Unknown"
    clean = str(val).strip().upper()
    if clean in ("CR", "COMPLETE RESPONSE"):
        return "CR"
    if clean in ("PR", "PARTIAL RESPONSE"):
        return "PR"
    if clean in ("SD", "STABLE DISEASE"):
        return "SD"
    if clean in ("PD", "PROGRESSIVE DISEASE", "PROGRESSION"):
        return "PD"
    return "Unknown"


def standardize_timepoint(val: object) -> str:
    """Standardize biopsy timepoint relative to therapy."""
    if val is None:
        return "Unknown"
    clean = str(val).strip().lower()
    if any(k in clean for k in ("pre", "baseline", "naive", "screening")):
        return "Pre"
    if any(k in clean for k in ("post", "progression", "resistant", "relapse")):
        return "Post"
    if any(k in clean for k in ("on", "during", "cycle")):
        return "On-Treatment"
    return "Unknown"


def harmonize_obs_metadata(obs_df: pd.DataFrame | pl.DataFrame) -> Result[pl.DataFrame, str]:
    """Harmonize AnnData .obs table into a unified Polars DataFrame."""
    try:
        df = pl.from_pandas(obs_df) if isinstance(obs_df, pd.DataFrame) else obs_df

        # Look for candidate response column
        cols = df.columns
        resp_col = next((c for c in cols if any(k in c.lower() for k in ("response", "recist", "clinical_benefit", "dcb"))), None)

        if resp_col:
            df = df.with_columns(
                pl.col(resp_col).map_elements(binarize_response, return_dtype=pl.Float64).alias("response_binary"),
                pl.col(resp_col).map_elements(standardize_recist, return_dtype=pl.String).alias("response_recist"),
            )

        # Look for candidate timepoint column
        tp_col = next((c for c in cols if any(k in c.lower() for k in ("timepoint", "time_point", "biopsy", "treatment_status"))), None)
        if tp_col:
            df = df.with_columns(
                pl.col(tp_col).map_elements(standardize_timepoint, return_dtype=pl.String).alias("biopsy_timepoint")
            )

        return Success(df)
    except Exception as exc:
        return Failure(f"Failed to harmonize metadata: {exc}")
