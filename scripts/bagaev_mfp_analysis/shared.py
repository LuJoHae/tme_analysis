"""Shared pure evaluation and table-I/O helpers for the Bagaev et al. (2021) MFP pipeline."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import TypeVar
import numpy as np
import polars as pl
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success

from tme_response import evaluate_subtype_response_associations, fit_univariable_cox

T = TypeVar("T")

SURVIVAL_ENDPOINTS: tuple[tuple[str, str, str], ...] = (
    ("OS", "os_months", "os_event"),
    ("PFS", "pfs_months", "pfs_event"),
)


def unwrap_or_raise(res: Result[T, str], context: str = "") -> T:
    """Unwrap a Result at an imperative CLI boundary or raise RuntimeError."""
    if isinstance(res, Failure):
        prefix = f"{context}: " if context else ""
        raise RuntimeError(f"{prefix}{res.failure()}")
    return res.unwrap()


def assert_inside_data_root(path: Path, data_root: Path) -> Path:
    """Ensure a resolved data path resides strictly inside `data_root` (e.g. `/storage/halu/data-test`)."""
    res_p, res_root = path.resolve(), data_root.resolve()
    if not res_p.is_relative_to(res_root):
        raise ValueError(f"Path '{res_p}' escapes configured data_root '{res_root}'.")
    return res_p


def write_table(df: pl.DataFrame, out_dir: Path, stem: str, *, write_tsv: bool = True) -> None:
    """Write a Polars DataFrame to `<stem>.parquet` and optionally `<stem>.tsv`."""
    df.write_parquet(out_dir / f"{stem}.parquet")
    if write_tsv:
        df.write_csv(out_dir / f"{stem}.tsv", separator="\t")


def eval_tagged_response(
    df: pl.DataFrame,
    *,
    subtype_col: str,
    cohort_name: str,
    sample_scope: str,
    method: str,
    stratify_by_cohort: bool,
    extra_exprs: Sequence[pl.Expr] = (),
) -> list[pl.DataFrame]:
    """Run `evaluate_subtype_response_associations` and return a 0-or-1-element tagged list."""
    res = evaluate_subtype_response_associations(
        df,
        subtype_col=subtype_col,
        cohort_col=Some("cohort") if stratify_by_cohort else Nothing,
    )
    if isinstance(res, Failure):
        return []
    return [
        res.unwrap().with_columns(
            [
                pl.lit(cohort_name).alias("cohort"),
                pl.lit(sample_scope).alias("sample_scope"),
                pl.lit(method).alias("method"),
                *extra_exprs,
            ]
        )
    ]


def evaluate_cox_contrasts(
    pre_df: pl.DataFrame,
    cohort_ids: Sequence[str],
    contrasts: Sequence[tuple[str, pl.Expr]],
) -> list[dict[str, str | float | int]]:
    """Evaluate unscaled binary OS/PFS Cox PH contrasts (cohort-stratified for PanCancer-iAtlas)."""
    scopes = [("PanCancer-iAtlas", pre_df), *((cid, pre_df.filter(pl.col("cohort") == cid)) for cid in cohort_ids)]
    rows: list[dict[str, str | float | int]] = []
    for endpoint, t_col, e_col in SURVIVAL_ENDPOINTS:
        for scope_id, sub_s in scopes:
            valid_s = sub_s.filter(
                pl.col(t_col).is_not_null()
                & pl.col(t_col).is_finite()
                & (pl.col(t_col) > 0.0)
                & pl.col(e_col).cast(pl.Float64).is_in([0.0, 1.0])
            )
            if valid_s.height < 15:
                continue
            times = valid_s[t_col].to_numpy().astype(np.float64).tolist()
            events = valid_s[e_col].to_numpy().astype(np.int64).tolist()
            if sum(events) < 5:
                continue
            strata_arg: Maybe[Sequence[str | int]] = (
                Some(valid_s["cohort"].cast(pl.String).to_list()) if scope_id == "PanCancer-iAtlas" else Nothing
            )
            for contrast_name, expr in contrasts:
                ind = valid_s.select(expr.cast(pl.Float64).alias("ind"))["ind"].to_list()
                cox_res = fit_univariable_cox(times, events, ind, standardize_score=False, strata=strata_arg)
                if isinstance(cox_res, Success):
                    cx = cox_res.unwrap()
                    rows.append(
                        {
                            "cohort": scope_id,
                            "endpoint": endpoint,
                            "contrast": contrast_name,
                            "n_samples": int(cx.n_samples),
                            "n_events": int(cx.n_events),
                            "hazard_ratio": float(cx.hazard_ratio),
                            "hr_ci95_low": float(cx.hr_ci_lower),
                            "hr_ci95_high": float(cx.hr_ci_upper),
                            "wald_p": float(cx.p_value),
                        }
                    )
    return rows
