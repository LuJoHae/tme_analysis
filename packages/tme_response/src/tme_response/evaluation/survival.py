"""Time-to-event survival evaluation (Harrell's C-index & Cox Proportional Hazards)."""

from __future__ import annotations

from typing import Sequence
import numpy as np
import pandas as pd
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from statsmodels.duration.hazard_regression import PHReg  # type: ignore[import-untyped]


class CIndexResult(BaseModel):
    """Concordance index evaluation result with bootstrap confidence intervals."""
    model_config = ConfigDict(frozen=True)

    c_index: float
    ci_lower: float
    ci_upper: float
    p_value_vs_half: float
    n_samples: int
    n_events: int
    n_comparable_pairs: int


class CoxHazardResult(BaseModel):
    """Univariable Cox Proportional Hazards regression result."""
    model_config = ConfigDict(frozen=True)

    hazard_ratio: float
    hr_ci_lower: float
    hr_ci_upper: float
    coefficient: float
    se_coefficient: float
    p_value: float
    n_samples: int
    n_events: int


def compute_c_index_raw(
    times: np.ndarray,
    events: np.ndarray,
    scores: np.ndarray,
    higher_is_better: bool = True,
) -> tuple[float, int]:
    """Calculate Harrell's C-index across all comparable pairs using vectorized upper-triangle broadcasting."""
    n = len(times)
    if n < 2:
        return 0.5, 0

    t_i = times[:, None]
    t_j = times[None, :]
    e_i = events[:, None]
    e_j = events[None, :]
    s_i = scores[:, None]
    s_j = scores[None, :]

    tri = np.triu(np.ones((n, n), dtype=bool), k=1)

    case1 = (t_i < t_j) & (e_i == 1) & tri
    case2 = (t_i > t_j) & (e_j == 1) & tri
    case3 = (t_i == t_j) & (e_i == 1) & (e_j == 1) & tri

    comparable = case1 | case2 | case3
    n_comp = int(np.sum(comparable))
    if n_comp == 0:
        return 0.5, 0

    diff = (s_j - s_i) if higher_is_better else (s_i - s_j)

    conc = float(
        np.sum(case1 & (diff > 0))
        + np.sum(case2 & (diff < 0))
        + np.sum(case3 & (diff == 0))
    )
    tied = float(
        np.sum(case1 & (diff == 0))
        + np.sum(case2 & (diff == 0))
        + np.sum(case3 & (diff != 0))
    )

    c = (conc + 0.5 * tied) / n_comp
    return float(c), n_comp


def calculate_c_index_bootstrap(
    times: Sequence[float],
    events: Sequence[float | int | bool],
    scores: Sequence[float],
    n_bootstrap: int = 500,
    alpha: float = 0.05,
    seed: int = 42,
    higher_is_better: bool = True,
) -> Result[CIndexResult, str]:
    """Calculate Harrell's C-index with empirical bootstrap 95% CI and two-sided p-value vs chance (0.50)."""
    try:
        t_arr = np.asarray(times, dtype=float)
        e_arr = np.asarray(events, dtype=int)
        s_arr = np.asarray(scores, dtype=float)

        valid_mask = np.isfinite(t_arr) & np.isfinite(e_arr) & np.isfinite(s_arr) & (t_arr > 0)
        t_clean = t_arr[valid_mask]
        e_clean = e_arr[valid_mask]
        s_clean = s_arr[valid_mask]

        n_samples = len(t_clean)
        n_events = int(np.sum(e_clean))

        if n_samples < 5 or n_events < 2:
            return Failure(f"Insufficient survival data: n_samples={n_samples}, n_events={n_events}")

        point_c, n_comp = compute_c_index_raw(t_clean, e_clean, s_clean, higher_is_better=higher_is_better)
        if n_comp == 0:
            return Failure("Zero comparable pairs in survival dataset.")

        # Bootstrap resampling
        rng = np.random.default_rng(seed)
        boot_c_list: list[float] = []

        for _ in range(n_bootstrap):
            idx = rng.choice(n_samples, size=n_samples, replace=True)
            t_b, e_b, s_b = t_clean[idx], e_clean[idx], s_clean[idx]
            if np.sum(e_b) < 1:
                continue
            b_c, b_comp = compute_c_index_raw(t_b, e_b, s_b, higher_is_better=higher_is_better)
            if b_comp > 0:
                boot_c_list.append(b_c)

        if len(boot_c_list) < 20:
            return Success(
                CIndexResult(
                    c_index=point_c,
                    ci_lower=point_c,
                    ci_upper=point_c,
                    p_value_vs_half=1.0,
                    n_samples=n_samples,
                    n_events=n_events,
                    n_comparable_pairs=n_comp,
                )
            )

        boot_arr = np.asarray(boot_c_list)
        ci_lower = float(np.percentile(boot_arr, 100.0 * (alpha / 2.0)))
        ci_upper = float(np.percentile(boot_arr, 100.0 * (1.0 - alpha / 2.0)))

        # Two-sided empirical p-value vs 0.50
        p_le = float(np.mean(boot_arr <= 0.5))
        p_ge = float(np.mean(boot_arr >= 0.5))
        p_val = min(1.0, 2.0 * min(p_le, p_ge))

        return Success(
            CIndexResult(
                c_index=point_c,
                ci_lower=ci_lower,
                ci_upper=ci_upper,
                p_value_vs_half=p_val,
                n_samples=n_samples,
                n_events=n_events,
                n_comparable_pairs=n_comp,
            )
        )
    except Exception as exc:
        return Failure(f"Failed to calculate C-index: {exc}")


def fit_univariable_cox(
    times: Sequence[float],
    events: Sequence[float | int | bool],
    scores: Sequence[float],
    standardize_score: bool = True,
    strata: Maybe[Sequence[str | int]] = Nothing,
) -> Result[CoxHazardResult, str]:
    """Fit a univariable Cox Proportional Hazards regression model using statsmodels PHReg.

    If standardize_score is True, scores are z-scored so the Hazard Ratio represents
    the hazard multiplier per 1 standard deviation increase in score.
    When `strata` is Some(...), fits a stratified Cox model allowing cohort-specific baseline hazards.
    """
    try:
        t_arr = np.asarray(times, dtype=float)
        e_arr = np.asarray(events, dtype=int)
        s_arr = np.asarray(scores, dtype=float)

        valid_mask = np.isfinite(t_arr) & np.isfinite(e_arr) & np.isfinite(s_arr) & (t_arr > 0)
        t_clean = t_arr[valid_mask]
        e_clean = e_arr[valid_mask]
        s_clean = s_arr[valid_mask]

        n_samples = len(t_clean)
        n_events = int(np.sum(e_clean))

        if n_samples < 10 or n_events < 3:
            return Failure(f"Insufficient survival data for Cox model: n_samples={n_samples}, n_events={n_events}")

        if standardize_score:
            std_val = float(np.std(s_clean))
            if std_val > 1e-8:
                s_model = (s_clean - float(np.mean(s_clean))) / std_val
            else:
                return Failure("Score has zero variance; cannot fit Cox regression.")
        else:
            if float(np.std(s_clean)) <= 1e-8:
                return Failure("Score has zero variance; cannot fit Cox regression.")
            s_model = s_clean

        u_vals = np.unique(s_model)
        if len(u_vals) <= 2:
            if len(u_vals) < 2:
                return Failure("Binary predictor has fewer than 2 levels; cannot fit Cox regression.")
            arm0_mask = s_model == u_vals[0]
            arm1_mask = s_model == u_vals[1]
            if int(np.sum(arm0_mask)) < 2 or int(np.sum(arm1_mask)) < 2:
                return Failure("Binary Cox regression requires at least 2 samples per exposure arm.")
            if int(np.sum(e_clean[arm0_mask])) == 0 or int(np.sum(e_clean[arm1_mask])) == 0:
                return Failure("Zero-event arm in binary Cox regression indicates complete separation.")

        strata_arg: list[str] | None = None
        match strata:
            case Some(s_vals):
                if len(s_vals) != len(t_arr):
                    return Failure(f"Strata length ({len(s_vals)}) does not match times ({len(t_arr)}).")
                strata_arg = [str(v) for v, keep in zip(s_vals, valid_mask.tolist(), strict=True) if keep]
            case _:
                pass

        df_model = pd.DataFrame({"time": t_clean, "event": e_clean, "score": s_model})
        mod = (
            PHReg(
                endog=df_model["time"],
                exog=df_model[["score"]],
                status=df_model["event"],
                strata=np.asarray(strata_arg, dtype=object),
            )
            if strata_arg is not None
            else PHReg(endog=df_model["time"], exog=df_model[["score"]], status=df_model["event"])
        )
        res = mod.fit(disp=False)

        coef = float(res.params[0])
        bse = float(res.bse[0])
        p_val = float(res.pvalues[0])
        ci = res.conf_int()
        hr = float(np.exp(np.clip(coef, -50.0, 50.0)))
        hr_lower = float(np.exp(ci[0, 0]))
        hr_upper = float(np.exp(ci[0, 1]))

        if (
            not np.isfinite(coef)
            or not np.isfinite(bse)
            or not np.isfinite(p_val)
            or not np.isfinite(hr_lower)
            or not np.isfinite(hr_upper)
            or abs(coef) > 10.0
            or bse > 10.0
            or hr_upper > 1e4
        ):
            return Failure(
                f"Degenerate Cox PH fit (quasi-complete separation): coef={coef:.3f}, se={bse:.3f}."
            )

        return Success(
            CoxHazardResult(
                hazard_ratio=hr,
                hr_ci_lower=hr_lower,
                hr_ci_upper=hr_upper,
                coefficient=coef,
                se_coefficient=bse,
                p_value=p_val,
                n_samples=n_samples,
                n_events=n_events,
            )
        )
    except Exception as exc:
        return Failure(f"Failed to fit Cox Proportional Hazards model: {exc}")
