"""Cryptographic checksum calculation and file verification."""

from __future__ import annotations

import hashlib
from pathlib import Path
from returns.result import Failure, Result, Success

from ..models import ChecksumSpec


def compute_file_hash(
    path: Path,
    algorithm: str = "sha256",
    chunk_size: int = 1024 * 1024,
) -> Result[str, str]:
    """Compute cryptographic hash of a file using buffered chunked streaming."""
    return (
        Failure(f"File does not exist: {path}")
        if not path.is_file()
        else _stream_hash(path, algorithm, chunk_size)
    )


def _stream_hash(path: Path, algorithm: str, chunk_size: int) -> Result[str, str]:
    try:
        hasher = hashlib.new(algorithm)
        with open(path, "rb") as fh:
            while chunk := fh.read(chunk_size):
                hasher.update(chunk)
        return Success(hasher.hexdigest())
    except Exception as exc:
        return Failure(f"Failed to calculate {algorithm} for {path}: {exc}")


def verify_checksum(path: Path, spec: ChecksumSpec) -> Result[bool, str]:
    """Verify that a given file matches the expected cryptographic checksum specification."""
    match compute_file_hash(path, spec.algorithm):
        case Failure(err):
            return Failure(err)
        case Success(computed_hash):
            return (
                Success(True)
                if computed_hash.lower() == spec.expected_hash.lower()
                else Failure(
                    f"Checksum mismatch for {path.name}: "
                    f"expected {spec.expected_hash.lower()}, computed {computed_hash.lower()}"
                )
            )
