"""InstaPrism: Fast Bayesian deconvolution."""

from __future__ import annotations

__version__ = "0.1.0"

from instaprism._instaprism import (
    bayes_prism,
    insta_prism,
)

__all__ = [
    "__version__",
    "bayes_prism",
    "insta_prism",
]
