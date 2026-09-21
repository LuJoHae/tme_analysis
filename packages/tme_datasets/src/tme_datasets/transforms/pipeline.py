"""PyTorch-style monadic composition of AnnData transformations."""

from __future__ import annotations

from typing import Callable, Sequence
import anndata as ad
from returns.result import Result, Success

from ..types import PerturbationTransform

TransformFn = Callable[[ad.AnnData], Result[ad.AnnData, str]]


class ComposeTransforms:
    """Sequential pipeline of AnnData perturbation transforms.

    Evaluates transforms sequentially using monadic bind; short-circuits on any Failure.
    """

    def __init__(self, transforms: Sequence[TransformFn]) -> None:
        self.transforms = tuple(transforms)

    def __call__(self, adata: ad.AnnData) -> Result[ad.AnnData, str]:
        current_res: Result[ad.AnnData, str] = Success(adata)
        for t in self.transforms:
            current_res = current_res.bind(t)
        return current_res
