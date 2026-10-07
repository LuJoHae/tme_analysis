"""Scikit-learn compatible adapter for Stability Selection.

Allows StabilitySelector to be dropped into standard scikit-learn Pipelines,
GridSearchCV, and FeatureUnion.
"""

from typing import Sequence
import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin  # type: ignore[import-untyped]
from returns.maybe import Maybe, Some, Nothing
from returns.result import Success, Failure

from selective_inference.stability_selection.types import (
    StabilityParameters,
    StabilityResult,
    SamplingType,
    Assumption,
)
from selective_inference.stability_selection.bounds import resolve_stability_parameters
from selective_inference.stability_selection.core import (
    run_stability_selection,
    default_lasso_path_fitter,
)


class StabilitySelector(BaseEstimator, TransformerMixin):  # type: ignore[misc]
    """Scikit-learn Transformer for Stability Selection.

    Parameters
    ----------
    cutoff : float, optional
        Threshold probability pi_thr (e.g., 0.75).
    q : float, optional
        Expected number of selected variables per subsample.
    pfer : float, optional
        Per-Family Error Rate (expected false positive discoveries).
    B : int, default 50
        Number of complementary pairs (or subsamples). Total subsamples = 2 * B.
    weakness : float, default 0.5
        Randomized Lasso feature weighting parameter alpha in (0, 1].
    sampling_type : {"SS", "MB"}, default "SS"
        "SS": Complementary pairs (Shah & Samworth 2013).
        "MB": Subsampling without replacement (Meinshausen & Bühlmann 2010).
    assumption : {"unimodal", "none"}, default "unimodal"
        Error bound assumption.
    lambdas : Sequence[float], optional
        Regularization path penalties.
    random_state : int, default 42
        Seed for reproducibility.
    """

    def __init__(
        self,
        cutoff: float | None = 0.75,
        q: float | None = None,
        pfer: float | None = 1.0,
        B: int = 50,
        weakness: float = 0.5,
        sampling_type: SamplingType = "SS",
        assumption: Assumption = "unimodal",
        lambdas: Sequence[float] | None = None,
        random_state: int = 42,
        n_jobs: int = 1,
    ) -> None:
        self.cutoff = cutoff
        self.q = q
        self.pfer = pfer
        self.B = B
        self.weakness = weakness
        self.sampling_type = sampling_type
        self.assumption = assumption
        self.lambdas = lambdas
        self.random_state = random_state
        self.n_jobs = n_jobs

    def fit(self, X: np.ndarray, y: np.ndarray) -> "StabilitySelector":
        """Fit Stability Selection on data X and target y."""
        X_arr = np.asarray(X, dtype=np.float64)
        y_arr = np.asarray(y, dtype=np.float64)

        n_samples, n_features = X_arr.shape

        cutoff_m = Some(self.cutoff) if self.cutoff is not None else Nothing
        q_m = Some(self.q) if self.q is not None else Nothing
        pfer_m = Some(self.pfer) if self.pfer is not None else Nothing

        # Resolve stability parameters
        param_res = resolve_stability_parameters(
            p=n_features,
            cutoff=cutoff_m,
            q=q_m,
            pfer=pfer_m,
            B=self.B,
            sampling_type=self.sampling_type,
            assumption=self.assumption,
        )

        match param_res:
            case Failure(err):
                raise ValueError(f"Stability parameter resolution failed: {err}")
            case Success(params):
                resolved_params = params

        # Run stability selection
        result_res = run_stability_selection(
            X=X_arr,
            y=y_arr,
            parameters=resolved_params,
            lambdas=self.lambdas,
            fitter=default_lasso_path_fitter,
            weakness=self.weakness,
            seed=self.random_state,
            n_jobs=self.n_jobs,
        )

        match result_res:
            case Failure(err):
                raise RuntimeError(f"Stability selection failed: {err}")
            case Success(res):
                self.result_ = res
                self.stability_scores_ = np.asarray(res.stability_scores, dtype=np.float64)
                self.stability_matrix_ = np.asarray(res.stability_matrix, dtype=np.float64)
                self.selected_indices_ = np.asarray(res.selected_indices, dtype=np.int64)
                self.n_features_in_ = n_features
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        """Select features that met the stability cutoff threshold."""
        if not hasattr(self, "selected_indices_"):
            raise ValueError("StabilitySelector must be fitted before calling transform.")
        X_arr = np.asarray(X)
        return X_arr[:, self.selected_indices_]

    def get_support(self, indices: bool = False) -> np.ndarray:
        """Get a mask, or integer index list, of the features selected."""
        if not hasattr(self, "selected_indices_"):
            raise ValueError("StabilitySelector must be fitted before calling get_support.")
        if indices:
            return self.selected_indices_
        mask = np.zeros(self.n_features_in_, dtype=bool)
        mask[self.selected_indices_] = True
        return mask
