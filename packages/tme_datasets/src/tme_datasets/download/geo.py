"""Zero-dependency automated downloaders for NCBI GEO and EMBL-EBI ArrayExpress datasets."""

from __future__ import annotations

from ftplib import FTP
from pathlib import Path
from typing import Sequence
import urllib.error
import urllib.request
from returns.result import Failure, Result, Success

from ..logging import get_logger
from .fetcher import download_single_file

logger = get_logger("download.geo")

NCBI_GEO_FTP_HOST = "ftp.ncbi.nlm.nih.gov"
EBI_FTP_HOST = "ftp.ebi.ac.uk"


def get_geo_suppl_url(gse_id: str, filename: str) -> str:
    """Build canonical NCBI HTTPS URL for a supplementary file."""
    prefix = gse_id[:-3]
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}nnn/{gse_id}/suppl/{filename}"


def list_geo_supplementary_files(gse_id: str) -> Result[tuple[str, ...], str]:
    """List supplementary files available for a GSE accession via anonymous FTP."""
    prefix = gse_id[:-3]
    remote_dir = f"/geo/series/{prefix}nnn/{gse_id}/suppl"
    try:
        logger.debug("Connecting to NCBI FTP (%s) to list %s...", NCBI_GEO_FTP_HOST, remote_dir)
        ftp = FTP(NCBI_GEO_FTP_HOST, timeout=30)
        ftp.login()
        ftp.cwd(remote_dir)
        files = tuple(ftp.nlst())
        ftp.quit()
        return Success(files)
    except Exception as exc:
        msg = f"Failed to list GEO supplementary files for {gse_id} via FTP: {exc}"
        logger.warning(msg)
        return Failure(msg)


def download_geo_supplementary(
    gse_id: str,
    dest_dir: Path,
    expected_files: Sequence[str] | None = None,
) -> Result[tuple[Path, ...], str]:
    """Download supplementary files for a given NCBI GEO GSE ID.

    Tries HTTPS download for specified or listed files first, with FTP fallback.
    """
    dest_dir.mkdir(parents=True, exist_ok=True)

    # Determine files to download
    target_files: tuple[str, ...]
    if expected_files is not None and len(expected_files) > 0:
        target_files = tuple(expected_files)
    else:
        list_res = list_geo_supplementary_files(gse_id)
        match list_res:
            case Success(files):
                target_files = files
            case Failure(err):
                # If FTP listing fails, check if files are already on disk
                existing = [p for p in dest_dir.iterdir() if p.is_file() and p.stat().st_size > 0]
                if existing:
                    logger.info("Found %d existing files in %s, skipping download.", len(existing), dest_dir)
                    return Success(tuple(existing))
                return Failure(f"Cannot determine supplementary files for {gse_id}: {err}")

def download_single_geo_file(
    gse_id: str,
    fname: str,
    dest_dir: Path,
) -> Result[Path, str]:
    """Download a single GEO supplementary file with HTTPS and FTP fallback."""
    dest_file = dest_dir / fname
    if dest_file.exists() and dest_file.stat().st_size > 0:
        logger.debug("File already exists: %s", dest_file.name)
        return Success(dest_file)

    https_url = get_geo_suppl_url(gse_id, fname)
    match download_single_file(https_url, dest_file):
        case Success(p):
            return Success(p)
        case Failure(https_err):
            logger.warning("HTTPS download failed for %s (%s). Attempting FTP fallback...", fname, https_err)
            prefix = gse_id[:-3]
            try:
                ftp = FTP(NCBI_GEO_FTP_HOST, timeout=60)
                ftp.login()
                ftp.cwd(f"/geo/series/{prefix}nnn/{gse_id}/suppl")
                with open(dest_file, "wb") as f_out:
                    ftp.retrbinary(f"RETR {fname}", f_out.write)
                ftp.quit()
                return Success(dest_file)
            except Exception as ftp_err:
                msg = f"Failed to download {fname} via both HTTPS and FTP: {ftp_err}"
                logger.error(msg)
                return Failure(msg)


def download_geo_supplementary(
    gse_id: str,
    dest_dir: Path,
    expected_files: Sequence[str] | None = None,
    max_workers: int = 4,
) -> Result[tuple[Path, ...], str]:
    """Download supplementary files for a given NCBI GEO GSE ID in parallel.

    Tries HTTPS download for specified or listed files first, with FTP fallback.
    """
    import concurrent.futures

    dest_dir.mkdir(parents=True, exist_ok=True)

    # Determine files to download
    target_files: tuple[str, ...]
    if expected_files is not None and len(expected_files) > 0:
        target_files = tuple(expected_files)
    else:
        list_res = list_geo_supplementary_files(gse_id)
        match list_res:
            case Success(files):
                target_files = files
            case Failure(err):
                existing = [p for p in dest_dir.iterdir() if p.is_file() and p.stat().st_size > 0]
                if existing:
                    logger.info("Found %d existing files in %s, skipping download.", len(existing), dest_dir)
                    return Success(tuple(existing))
                return Failure(f"Cannot determine supplementary files for {gse_id}: {err}")

    workers = min(max_workers, len(target_files)) if target_files else 1
    if workers > 1:
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            future_to_file = {
                executor.submit(download_single_geo_file, gse_id, fname, dest_dir): fname
                for fname in target_files
            }
            downloaded: list[Path] = []
            for future in concurrent.futures.as_completed(future_to_file):
                match future.result():
                    case Success(p):
                        downloaded.append(p)
                    case Failure(err):
                        return Failure(err)
            return Success(tuple(downloaded))

    downloaded_seq: list[Path] = []
    for fname in target_files:
        match download_single_geo_file(gse_id, fname, dest_dir):
            case Success(p):
                downloaded_seq.append(p)
            case Failure(err):
                return Failure(err)
    return Success(tuple(downloaded_seq))


def download_arrayexpress_files(
    accession: str,
    dest_dir: Path,
    expected_files: Sequence[str] | None = None,
) -> Result[tuple[Path, ...], str]:
    """Download files from an EMBL-EBI ArrayExpress / BioStudies dataset."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    suffix = accession[-3:]
    prefix_dir = f"{accession.split('-')[0]}-{accession.split('-')[1]}-" if "-" in accession else "E-MTAB-"
    remote_ftp_dir = f"/biostudies/fire/{prefix_dir}/{suffix}/{accession}/Files"

    target_files: tuple[str, ...]
    if expected_files is not None and len(expected_files) > 0:
        target_files = tuple(expected_files)
    else:
        try:
            logger.debug("Connecting to EBI FTP (%s) to list %s...", EBI_FTP_HOST, remote_ftp_dir)
            ftp = FTP(EBI_FTP_HOST, timeout=30)
            ftp.login()
            ftp.cwd(remote_ftp_dir)
            target_files = tuple(ftp.nlst())
            ftp.quit()
        except Exception as exc:
            existing = [p for p in dest_dir.iterdir() if p.is_file() and p.stat().st_size > 0]
            if existing:
                return Success(tuple(existing))
            return Failure(f"Cannot list ArrayExpress files for {accession}: {exc}")

    downloaded: list[Path] = []
    for fname in target_files:
        dest_file = dest_dir / fname
        if dest_file.exists() and dest_file.stat().st_size > 0:
            downloaded.append(dest_file)
            continue

        # Try EBI BioStudies HTTPS URL
        https_url = f"https://www.ebi.ac.uk/biostudies/files/{accession}/{fname}"
        match download_single_file(https_url, dest_file):
            case Success(p):
                downloaded.append(p)
            case Failure(https_err):
                logger.warning("HTTPS download failed for %s. Attempting EBI FTP fallback...", fname)
                try:
                    ftp = FTP(EBI_FTP_HOST, timeout=60)
                    ftp.login()
                    ftp.cwd(remote_ftp_dir)
                    with open(dest_file, "wb") as f_out:
                        ftp.retrbinary(f"RETR {fname}", f_out.write)
                    ftp.quit()
                    downloaded.append(dest_file)
                except Exception as ftp_err:
                    msg = f"Failed to download {fname} from ArrayExpress via both HTTPS and FTP: {ftp_err}"
                    logger.error(msg)
                    return Failure(msg)

    return Success(tuple(downloaded))
