"""Adapters module for third-party deconvolution frameworks."""

from bayesprism.adapters.rectangle import (
    RectangleConfig,
    RectangleDeconvResult,
    create_signature_from_matrix,
    build_rectangle_signatures,
    deconvolve_rectangle,
    is_rectangle_available,
)

from bayesprism.adapters.cibersort import (
    CibersortConfig,
    CibersortDeconvResult,
    deconvolve_cibersort,
)

from bayesprism.adapters.cibersortx_docker import (
    CibersortXDockerConfig,
    CibersortXResult,
    run_cibersortx_docker,
    is_docker_available,
    format_cibersort_matrix,
    parse_cibersortx_results,
    build_cibersortx_docker_command,
    resolve_cibersortx_credentials,
)

__all__ = [
    "RectangleConfig",
    "RectangleDeconvResult",
    "create_signature_from_matrix",
    "build_rectangle_signatures",
    "deconvolve_rectangle",
    "is_rectangle_available",
    "CibersortConfig",
    "CibersortDeconvResult",
    "deconvolve_cibersort",
    "CibersortXDockerConfig",
    "CibersortXResult",
    "run_cibersortx_docker",
    "is_docker_available",
    "format_cibersort_matrix",
    "parse_cibersortx_results",
    "build_cibersortx_docker_command",
    "resolve_cibersortx_credentials",
]

