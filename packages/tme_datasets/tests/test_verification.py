"""Unit tests for cryptographic verification."""

from pathlib import Path
from returns.result import Success
from tme_datasets.download.verification import compute_file_hash, verify_checksum
from tme_datasets.models import ChecksumSpec


def test_compute_file_hash(tmp_path: Path) -> None:
    test_file = tmp_path / "test.txt"
    test_file.write_text("TME analysis test content", encoding="utf-8")

    res = compute_file_hash(test_file, algorithm="sha256")
    assert isinstance(res, Success)
    hash_str = res.unwrap()
    assert len(hash_str) == 64

    # Verification success
    spec = ChecksumSpec(algorithm="sha256", expected_hash=hash_str)
    v_res = verify_checksum(test_file, spec)
    assert isinstance(v_res, Success)
    assert v_res.unwrap() is True


def test_verify_checksum_mismatch(tmp_path: Path) -> None:
    test_file = tmp_path / "test.txt"
    test_file.write_text("Hello World", encoding="utf-8")

    spec = ChecksumSpec(algorithm="sha256", expected_hash="0" * 64)
    v_res = verify_checksum(test_file, spec)
    assert not isinstance(v_res, Success)
