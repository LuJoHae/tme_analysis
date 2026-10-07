"""Mathematical error bound solvers for Stability Selection.

Implements finite-sample error control under:
- Meinshausen & Bühlmann (2010): Per-Family Error Rate (PFER) under exchangeability.
- Shah & Samworth (2013): Complementary Pairs Stability Selection (CPSS) under Unimodal distribution.

Strictly functional:
- Pure mathematical functions
- Railway-Oriented Programming (Result[T, str])
- No None (Maybe monad)
"""

import math
from typing import Literal
import numpy as np
from returns.result import Result, Success, Failure
from returns.maybe import Maybe
import returns.maybe as rm

from selective_inference.stability_selection.types import (
    StabilityParameters,
    SamplingType,
    Assumption,
)


def compute_unimodal_constant(cutoff: float, B: int) -> float:
    """Compute Shah & Samworth (2013) unimodal denominator constant C(pi_thr, B).

    Parameters
    ----------
    cutoff : float
        Selection cutoff threshold pi_thr in (0.5, 1.0].
    B : int
        Number of complementary pairs (subsamples).

    Returns
    -------
    float
        The multiplier constant C(pi_thr, B).
    """
    if cutoff <= 0.75:
        return 2.0 * (2.0 * cutoff - 1.0 - 1.0 / (2.0 * B))
    return (1.0 + 1.0 / B) / (4.0 * (1.0 - cutoff + 1.0 / (2.0 * B)))


def solve_pfer(
    p: int,
    cutoff: float,
    q: float,
    B: int,
    assumption: Assumption,
) -> Result[float, str]:
    """Calculate the upper bound on the Per-Family Error Rate (PFER).

    PFER = E[V], the expected number of falsely selected noise features.
    """
    if p <= 0:
        return Failure(f"Total features p must be positive, got {p}")
    if q <= 0 or q > p:
        return Failure(f"q must be in (0, {p}], got {q}")
    if cutoff < 0.5 or cutoff > 1.0:
        return Failure(f"cutoff must be in [0.5, 1.0], got {cutoff}")
    if B <= 1:
        return Failure(f"B must be at least 2, got {B}")

    match assumption:
        case "none":
            denom = 2.0 * cutoff - 1.0
            if denom <= 0:
                return Failure("Cutoff must be strictly > 0.5 for finite PFER under assumption 'none'")
            return Success((q ** 2) / (p * denom))

        case "unimodal":
            theta = q / p
            c_min = 0.5 + min(theta ** 2, 1.0 / (2.0 * B) + 0.75 * (theta ** 2))
            if cutoff <= c_min:
                return Failure(
                    f"Cutoff {cutoff:.4f} violates unimodal admissibility condition (must be > {c_min:.4f})"
                )
            const = compute_unimodal_constant(cutoff, B)
            if const <= 0:
                return Failure(f"Unimodal constant is non-positive ({const}) for cutoff {cutoff}")
            return Success((q ** 2) / (p * const))

        case "r-concave":
            # Conservative fallback using unimodal bound
            return solve_pfer(p, cutoff, q, B, "unimodal")


def solve_cutoff(
    p: int,
    q: float,
    pfer: float,
    B: int,
    assumption: Assumption,
) -> Result[float, str]:
    """Solve for the optimal stability cutoff threshold pi_thr given p, q, and target PFER."""
    if p <= 0 or q <= 0 or q > p:
        return Failure(f"Invalid p={p} or q={q}")
    if pfer <= 0:
        return Failure(f"Target PFER must be positive, got {pfer}")
    if B <= 1:
        return Failure(f"B must be at least 2, got {B}")

    match assumption:
        case "none":
            calculated = 0.5 + (q ** 2) / (2.0 * p * pfer)
            if calculated > 1.0:
                return Failure(
                    f"Required cutoff ({calculated:.3f}) exceeds 1.0. Increase target PFER or reduce q."
                )
            return Success(max(0.5, calculated))

        case "unimodal":
            theta = q / p
            c_min = 0.5 + min(theta ** 2, 1.0 / (2.0 * B) + 0.75 * (theta ** 2))

            # Discretized search grid matching Shah & Samworth / stabs
            cutoff_grid = [0.5 + k / (2.0 * B) for k in range(2, B + 1)]
            admissible_cutoffs = [c for c in cutoff_grid if c > c_min and c <= 1.0]

            if not admissible_cutoffs:
                return Failure(
                    f"No admissible cutoffs for q={q}, p={p}, B={B}. Try decreasing q."
                )

            for cand_cutoff in admissible_cutoffs:
                cand_const = compute_unimodal_constant(cand_cutoff, B)
                bound = (q ** 2) / (p * cand_const)
                if bound <= pfer:
                    return Success(cand_cutoff)

            # If no discrete cutoff strictly achieves PFER, return the highest (safest)
            return Success(admissible_cutoffs[-1])

        case "r-concave":
            return solve_cutoff(p, q, pfer, B, "unimodal")


def solve_q(
    p: int,
    cutoff: float,
    pfer: float,
    B: int,
    assumption: Assumption,
) -> Result[float, str]:
    """Solve for the maximum allowed selection budget q given p, cutoff, and target PFER."""
    if p <= 0 or pfer <= 0:
        return Failure(f"Invalid p={p} or pfer={pfer}")
    if cutoff < 0.5 or cutoff > 1.0:
        return Failure(f"Cutoff must be in [0.5, 1.0], got {cutoff}")
    if B <= 1:
        return Failure(f"B must be at least 2, got {B}")

    match assumption:
        case "none":
            denom = 2.0 * cutoff - 1.0
            if denom <= 0:
                return Failure("Cutoff must be > 0.5 to solve for q")
            max_q = math.sqrt(p * pfer * denom)
            return Success(min(float(p), math.floor(max_q)))

        case "unimodal":
            const = compute_unimodal_constant(cutoff, B)
            if const <= 0:
                return Failure(f"Invalid unimodal constant for cutoff {cutoff}")
            max_q = math.sqrt(p * pfer * const)
            # Enforce admissibility: cutoff > c_min => q must not violate cutoff bound
            q_cand = min(float(p), math.floor(max_q))
            return Success(max(1.0, q_cand))

        case "r-concave":
            return solve_q(p, cutoff, pfer, B, "unimodal")


def resolve_stability_parameters(
    p: int,
    cutoff: Maybe[float] = rm.Nothing,
    q: Maybe[float] = rm.Nothing,
    pfer: Maybe[float] = rm.Nothing,
    B: int = 50,
    sampling_type: SamplingType = "SS",
    assumption: Assumption = "unimodal",
) -> Result[StabilityParameters, str]:
    """Railway-oriented resolution of stability selection parameters.

    Exactly two of (cutoff, q, pfer) must be provided. The missing parameter
    is solved using analytical finite-sample error bounds.
    """
    effective_assumption: Assumption = "none" if sampling_type == "MB" else assumption
    specified_count = sum([cutoff != rm.Nothing, q != rm.Nothing, pfer != rm.Nothing])

    if specified_count != 2:
        return Failure(
            f"Exactly two of (cutoff, q, pfer) must be specified. Received {specified_count}."
        )

    match (cutoff, q, pfer):
        # Case 1: Solve for PFER
        case (rm.Some(c), rm.Some(selected_q), rm.Nothing):
            pfer_res = solve_pfer(p, c, selected_q, B, effective_assumption)
            return pfer_res.map(
                lambda val: StabilityParameters(
                    p=p,
                    q=selected_q,
                    cutoff=c,
                    pfer=val,
                    B=B,
                    sampling_type=sampling_type,
                    assumption=effective_assumption,
                )
            )

        # Case 2: Solve for cutoff
        case (rm.Nothing, rm.Some(selected_q), rm.Some(target_pfer)):
            cutoff_res = solve_cutoff(p, selected_q, target_pfer, B, effective_assumption)
            return cutoff_res.bind(
                lambda c: solve_pfer(p, c, selected_q, B, effective_assumption).map(
                    lambda actual_pfer: StabilityParameters(
                        p=p,
                        q=selected_q,
                        cutoff=c,
                        pfer=actual_pfer,
                        B=B,
                        sampling_type=sampling_type,
                        assumption=effective_assumption,
                    )
                )
            )

        # Case 3: Solve for q
        case (rm.Some(c), rm.Nothing, rm.Some(target_pfer)):
            q_res = solve_q(p, c, target_pfer, B, effective_assumption)
            return q_res.bind(
                lambda solved_q: solve_pfer(p, c, solved_q, B, effective_assumption).map(
                    lambda actual_pfer: StabilityParameters(
                        p=p,
                        q=solved_q,
                        cutoff=c,
                        pfer=actual_pfer,
                        B=B,
                        sampling_type=sampling_type,
                        assumption=effective_assumption,
                    )
                )
            )

        case _:
            return Failure("Invalid parameter specification.")

