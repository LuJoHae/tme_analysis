"""Automated URL downloading with resumable HTTP Range retries, checksums, and archive extraction."""

from __future__ import annotations

import http.client
from pathlib import Path
import tarfile
import time
import urllib.error
import urllib.request
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success

from ..logging import get_logger
from ..models import ChecksumSpec
from ..paths import is_cscs_environment
from .verification import verify_checksum

logger = get_logger("download")

CSCS_DOWNLOAD_BLOCKED_MSG = (
    "Internet dataset downloads are strictly disabled on CSCS Eiger. "
    "Download raw/preprocessed datasets on 'olm' first, transfer them manually to "
    "/capstor/store/cscs/empa/em25/halu/data-test via Globus, and run 'make stage-in'."
)


def _parse_total_bytes(
    status_code: int,
    offset: int,
    content_length_hdr: Maybe[str],
    content_range_hdr: Maybe[str],
) -> Maybe[int]:
    """Purely extract expected total file size in bytes from HTTP 200/206 headers."""
    match content_range_hdr:
        case Some(cr_val) if "/" in cr_val:
            total_part = cr_val.rsplit("/", 1)[-1].strip()
            if total_part.isdigit():
                return Some(int(total_part))
        case _:
            pass

    match content_length_hdr:
        case Some(cl_val) if cl_val.strip().isdigit():
            cl_int = int(cl_val.strip())
            return Some(offset + cl_int if status_code == 206 else cl_int)
        case _:
            return Nothing


def _verify_existing_destination(
    dest_path: Path,
    checksum: Maybe[ChecksumSpec],
) -> Result[bool, str]:
    """Verify an already-present destination file against optional size and hash specifications."""
    if not dest_path.is_file() or dest_path.stat().st_size == 0:
        return Success(False)

    actual_size = dest_path.stat().st_size
    match checksum:
        case Some(spec):
            match spec.size_bytes:
                case Some(expected_sz) if actual_size != expected_sz:
                    return Failure(
                        f"Existing file {dest_path.name} has size {actual_size}, expected {expected_sz}"
                    )
                case _:
                    pass
            return verify_checksum(dest_path, spec)
        case _:
            return Success(False)


def _stream_attempt(
    url: str,
    dest_path: Path,
    part_path: Path,
    resume: bool,
    timeout_sec: float,
) -> Result[tuple[int, Maybe[int], float], str]:
    """Execute a single HTTP GET (or Range-resumed GET) streaming attempt into part_path."""
    existing_offset = (
        part_path.stat().st_size if (resume and part_path.is_file()) else 0
    )
    headers: dict[str, str] = {"User-Agent": "Mozilla/5.0 (Python tme_datasets)"}
    if existing_offset > 0:
        headers["Range"] = f"bytes={existing_offset}-"

    req = urllib.request.Request(url, headers=headers)
    start_time = time.time()

    try:
        with urllib.request.urlopen(req, timeout=timeout_sec) as response:
            status_code = int(getattr(response, "status", 200))
            raw_cl = response.info().get("Content-Length")
            raw_cr = response.info().get("Content-Range")
            cl_hdr: Maybe[str] = Some(str(raw_cl)) if raw_cl is not None else Nothing
            cr_hdr: Maybe[str] = Some(str(raw_cr)) if raw_cr is not None else Nothing

            effective_offset = existing_offset if status_code == 206 else 0
            total_bytes_maybe = _parse_total_bytes(
                status_code=status_code,
                offset=effective_offset,
                content_length_hdr=cl_hdr,
                content_range_hdr=cr_hdr,
            )

            # Fast return if dest_path already matches remote Content-Length when no checksum was given
            if dest_path.is_file() and dest_path.stat().st_size > 0 and existing_offset == 0:
                match total_bytes_maybe:
                    case Some(remote_total) if dest_path.stat().st_size == remote_total:
                        logger.debug(
                            "File already present and size-verified: %s (%.1f MB)",
                            dest_path.name,
                            dest_path.stat().st_size / (1024 * 1024),
                        )
                        return Success((-1, total_bytes_maybe, 0.01))
                    case Some(remote_total):
                        logger.warning(
                            "Existing file %s is incomplete (%d / %d bytes). Re-downloading...",
                            dest_path.name,
                            dest_path.stat().st_size,
                            remote_total,
                        )
                        dest_path.unlink(missing_ok=True)
                    case _:
                        return Success((-1, Nothing, 0.01))

            open_mode = "ab" if (status_code == 206 and effective_offset > 0) else "wb"
            if effective_offset > 0:
                logger.info(
                    "Resuming %s from %.1f MB (HTTP %d) -> %s",
                    url,
                    effective_offset / (1024 * 1024),
                    status_code,
                    dest_path,
                )
            else:
                logger.info("Downloading %s -> %s", url, dest_path)

            chunk_size = 1024 * 1024  # 1 MiB
            downloaded = effective_offset
            last_reported_mb = downloaded / (1024 * 1024)

            with open(part_path, open_mode) as out_file:
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    out_file.write(chunk)
                    downloaded += len(chunk)
                    downloaded_mb = downloaded / (1024 * 1024)

                    if downloaded_mb - last_reported_mb >= 50.0:
                        match total_bytes_maybe:
                            case Some(tot_b) if tot_b > 0:
                                pct = (downloaded / tot_b) * 100.0
                                logger.info(
                                    "Downloading %s: %.1f / %.1f MB (%.1f%%)",
                                    dest_path.name,
                                    downloaded_mb,
                                    tot_b / (1024 * 1024),
                                    pct,
                                )
                            case _:
                                logger.info(
                                    "Downloading %s: %.1f MB downloaded",
                                    dest_path.name,
                                    downloaded_mb,
                                )
                        last_reported_mb = downloaded_mb

            elapsed = max(0.01, time.time() - start_time)
            return Success((downloaded, total_bytes_maybe, elapsed))
    except urllib.error.HTTPError as http_err:
        if http_err.code == 416 and part_path.exists():
            logger.warning(
                "HTTP 416 Range Not Satisfiable for %s; resetting partial file and retrying from byte 0.",
                part_path.name,
            )
            part_path.unlink(missing_ok=True)
        return Failure(f"HTTP {http_err.code}: {http_err.reason}")
    except (
        urllib.error.URLError,
        http.client.HTTPException,
        TimeoutError,
        ConnectionError,
        OSError,
    ) as net_err:
        return Failure(str(net_err))


def download_single_file(
    url: str,
    dest_path: Path,
    *,
    resume: bool = True,
    max_retries: int = 5,
    timeout_sec: float = 120.0,
    checksum: Maybe[ChecksumSpec] = Nothing,
) -> Result[Path, str]:
    """Download a file from a URL to dest_path with optional HTTP Range resumption and hash verification.

    Streams chunks to a temporary `.part` file, automatically resumes from partial bytes via HTTP 206
    upon transient network errors, optionally verifies size/MD5/SHA256, and atomically renames into place.
    """
    if is_cscs_environment():
        logger.error("%s (blocked URL: %s)", CSCS_DOWNLOAD_BLOCKED_MSG, url)
        return Failure(f"{CSCS_DOWNLOAD_BLOCKED_MSG} (blocked URL: {url})")

    try:
        dest_path.parent.mkdir(parents=True, exist_ok=True)

        # 1. If checksum is specified and dest_path already passes verification, return immediately
        match checksum:
            case Some(spec) if dest_path.is_file() and dest_path.stat().st_size > 0:
                match _verify_existing_destination(dest_path, checksum):
                    case Success(True):
                        logger.info(
                            "Verified existing file %s (%.1f MB, %s OK).",
                            dest_path.name,
                            dest_path.stat().st_size / (1024 * 1024),
                            spec.algorithm,
                        )
                        return Success(dest_path)
                    case Failure(v_err):
                        logger.warning("%s. Re-downloading...", v_err)
                        dest_path.unlink(missing_ok=True)
                    case _:
                        pass
            case _:
                pass

        part_path = dest_path.with_name(dest_path.name + ".part")
        if not resume and part_path.exists():
            part_path.unlink(missing_ok=True)

        # 2. If a complete `.part` file already matches expected size + checksum, promote it without network I/O
        if resume and part_path.is_file():
            match checksum:
                case Some(spec):
                    match spec.size_bytes:
                        case Some(exp_sz) if part_path.stat().st_size == exp_sz:
                            logger.info(
                                "Found complete partial file %s (%d bytes); verifying %s checksum...",
                                part_path.name,
                                exp_sz,
                                spec.algorithm,
                            )
                            match verify_checksum(part_path, spec):
                                case Success(_):
                                    part_path.replace(dest_path)
                                    logger.info("Promoted verified %s -> %s", part_path.name, dest_path.name)
                                    return Success(dest_path)
                                case Failure(part_err):
                                    logger.warning("%s. Resetting corrupt partial file...", part_err)
                                    part_path.unlink(missing_ok=True)
                        case _:
                            pass
                case _:
                    pass

        attempts = max(1, max_retries)
        last_error = "Unknown download error"

        for attempt_idx in range(attempts):
            attempt_res = _stream_attempt(
                url=url,
                dest_path=dest_path,
                part_path=part_path,
                resume=resume,
                timeout_sec=timeout_sec,
            )
            match attempt_res:
                case Success((-1, _, _)):
                    return Success(dest_path)
                case Success((downloaded, total_maybe, elapsed)):
                    size_complete = (
                        downloaded == total_maybe.unwrap()
                        if isinstance(total_maybe, Some)
                        else True
                    )
                    if not size_complete:
                        expected_b = total_maybe.value_or(0)
                        last_error = (
                            f"Download incomplete for {dest_path.name}: "
                            f"{downloaded}/{expected_b} bytes received"
                        )
                        logger.warning(
                            "%s (attempt %d/%d)",
                            last_error,
                            attempt_idx + 1,
                            attempts,
                        )
                        if not resume:
                            part_path.unlink(missing_ok=True)
                    else:
                        match checksum:
                            case Some(spec):
                                match spec.size_bytes:
                                    case Some(exp_sz) if downloaded != exp_sz:
                                        part_path.unlink(missing_ok=True)
                                        return Failure(
                                            f"Size mismatch for {dest_path.name}: "
                                            f"expected {exp_sz} bytes, received {downloaded}"
                                        )
                                    case _:
                                        pass
                                logger.info(
                                    "Verifying %s checksum for %s...",
                                    spec.algorithm,
                                    dest_path.name,
                                )
                                match verify_checksum(part_path, spec):
                                    case Failure(hash_err):
                                        part_path.unlink(missing_ok=True)
                                        return Failure(hash_err)
                                    case Success(_):
                                        pass
                            case _:
                                pass

                        part_path.replace(dest_path)
                        final_mb = downloaded / (1024 * 1024)
                        rate = final_mb / max(0.01, elapsed)
                        logger.info(
                            "Successfully downloaded %s (%.1f MB in %.1fs, %.1f MB/s)",
                            dest_path.name,
                            final_mb,
                            elapsed,
                            rate,
                        )
                        return Success(dest_path)
                case Failure(err_str):
                    last_error = err_str
                    logger.warning(
                        "Download attempt %d/%d failed for %s: %s",
                        attempt_idx + 1,
                        attempts,
                        dest_path.name,
                        err_str,
                    )
                    if not resume and part_path.exists():
                        part_path.unlink(missing_ok=True)

            if attempt_idx + 1 < attempts:
                backoff_sec = min(30.0, float(2**attempt_idx))
                time.sleep(backoff_sec)

        if not resume and part_path.exists():
            part_path.unlink(missing_ok=True)
        final_msg = f"Failed to download {url} -> {dest_path}: {last_error}"
        logger.error(final_msg)
        return Failure(final_msg)
    except Exception as exc:  # noqa: BLE001
        part_path = dest_path.with_name(dest_path.name + ".part")
        if not resume and part_path.exists():
            part_path.unlink(missing_ok=True)
        msg = f"Failed to download {url} -> {dest_path}: {exc}"
        logger.error(msg)
        return Failure(msg)


def unpack_tar(tar_path: Path, extract_dir: Path) -> Result[Path, str]:
    """Extract a tar archive to the specified directory."""
    if not tar_path.is_file():
        msg = f"Tar file not found: {tar_path}"
        logger.error(msg)
        return Failure(msg)

    try:
        extract_dir.mkdir(parents=True, exist_ok=True)
        logger.info("Extracting archive %s -> %s...", tar_path.name, extract_dir)
        start_time = time.time()

        with tarfile.open(tar_path, "r:*") as tar:
            if hasattr(tarfile, "data_filter"):
                tar.extractall(path=extract_dir, filter="data")
            else:
                tar.extractall(path=extract_dir)

        elapsed = time.time() - start_time
        logger.info(
            "Successfully extracted %s in %.1fs",
            tar_path.name,
            elapsed,
        )
        return Success(extract_dir)
    except Exception as exc:  # noqa: BLE001
        msg = f"Failed to unpack tar {tar_path}: {exc}"
        logger.error(msg)
        return Failure(msg)
