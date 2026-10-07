"""Pure functional clinical response label noise operator."""

from __future__ import annotations

import numpy as np
import polars as pl
from .models import LabelNoiseConfig


def apply_label_noise(
    clinical_df: pl.DataFrame,
    config: LabelNoiseConfig,
) -> pl.DataFrame:
    """Randomly flip binary response annotations with probability noise_rate.

    Simulates pseudo-progression, atypical response kinetics, and inter-radiologist RECIST discrepancy.
    """
    if config.noise_rate <= 0.0 or "response_binary" not in clinical_df.columns:
        return clinical_df

    rng = np.random.default_rng(config.seed)
    n = len(clinical_df)
    rate = float(min(1.0, max(0.0, config.noise_rate)))

    flip_mask = rng.uniform(0.0, 1.0, size=n) < rate
    responses = clinical_df["response_binary"].to_numpy().copy()

    # Flip 1.0 -> 0.0 and 0.0 -> 1.0 for selected indices
    responses[flip_mask] = 1.0 - responses[flip_mask]

    return clinical_df.with_columns(
        pl.Series("response_binary", responses)
    )
