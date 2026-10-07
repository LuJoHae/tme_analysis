"""Types and immutable data models for Stability Selection.

Follows strict functional programming principles:
- Immutable frozen models
- Typeclasses via typing.Protocol
- No None allowed (use Maybe monad)
"""

from typing import Protocol, Literal, Sequence
import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict, Field
from returns.maybe import Maybe, Some, Nothing


class PathFitter(Protocol):
    """Typeclass Protocol for high-dimensional path fitters.

    Maps data (X, y) and penalty sequence to a boolean matrix of shape
    (n_lambdas, n_features) indicating active feature selection.
    """

    def __call__(
        self,
        X: np.ndarray,
        y: np.ndarray,
        lambdas: Sequence[float] | np.ndarray,
        subsample_indices: np.ndarray | None = None,
    ) -> np.ndarray: ...


SamplingType = Literal["SS", "MB"]
Assumption = Literal["unimodal", "none", "r-concave"]


class StabilityParameters(BaseModel):
    """Parameters and theoretical bounds for stability selection."""

    model_config = ConfigDict(frozen=True)

    p: int = Field(gt=0, description="Total number of candidate features")
    q: float = Field(gt=0, description="Expected / average number of selected features per subsample")
    cutoff: float = Field(ge=0.5, le=1.0, description="Selection threshold pi_thr in [0.5, 1.0]")
    pfer: float = Field(ge=0.0, description="Upper bound on Per-Family Error Rate (expected false positives)")
    B: int = Field(gt=1, description="Number of subsample repetitions (or complementary pairs)")
    sampling_type: SamplingType = Field(description="Subsampling scheme: SS (Complementary Pairs) or MB")
    assumption: Assumption = Field(description="Distributional assumption for error bounds")


class StabilityResult(BaseModel):
    """Immutable result of a Stability Selection run."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    feature_names: tuple[str, ...]
    lambdas: tuple[float, ...]
    stability_matrix: tuple[tuple[float, ...], ...]  # shape: (n_features, n_lambdas)
    stability_scores: tuple[float, ...]              # max over valid lambdas per feature
    selected_features: tuple[str, ...]
    selected_indices: tuple[int, ...]
    empirical_q: float
    parameters: StabilityParameters
    expected_model_sizes: tuple[float, ...] = ()
    lambda_cutoff: Maybe[float] = Nothing
    unrestricted_stability_scores: tuple[float, ...] = ()

    def to_polars(self, include_unrestricted: bool = False) -> pl.DataFrame:
        """Convert stability scores to a sorted Polars DataFrame."""
        is_selected = tuple(
            feature in self.selected_features for feature in self.feature_names
        )
        data: dict[str, list[object]] = {
            "feature": list(self.feature_names),
            "stability_score": list(self.stability_scores),
            "selected": list(is_selected),
        }
        if include_unrestricted and self.unrestricted_stability_scores:
            data["unrestricted_score"] = list(self.unrestricted_stability_scores)
        return pl.DataFrame(data).sort("stability_score", descending=True)

    def to_path_polars(self) -> pl.DataFrame:
        """Convert full stability paths into long-form Polars DataFrame for Altair plotting."""
        cutoff_val = self.lambda_cutoff.value_or(0.0)
        has_cutoff = isinstance(self.lambda_cutoff, Some)

        records = [
            {
                "feature": self.feature_names[i],
                "lambda": float(self.lambdas[j]),
                "log_lambda": float(np.log10(self.lambdas[j])) if self.lambdas[j] > 0 else 0.0,
                "selection_probability": float(self.stability_matrix[i][j]),
                "selected": self.feature_names[i] in self.selected_features,
                "expected_model_size": float(self.expected_model_sizes[j])
                if j < len(self.expected_model_sizes)
                else 0.0,
                "in_budget": bool(self.lambdas[j] >= cutoff_val) if has_cutoff else True,
            }
            for i in range(len(self.feature_names))
            for j in range(len(self.lambdas))
        ]
        return pl.DataFrame(records)

