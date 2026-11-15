"""GTDB-Tk 2.7.1 R232 reference cache lifecycle."""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

from .cache_runtime import (
    CacheIntegrityError,
    append_log,
    ensure_fresh_directory_cache,
    readable_nonempty_file,
    resolved_path,
    run_logged,
)
from .download_retries import retry_download

GTDBTK_RELEASE = "r232"
GTDBTK_DATABASE_URL = (
    "https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/"
    "gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz"
)
GTDBTK_REQUIRED_FILES = (
    "metadata/metadata.txt",
    "taxonomy/gtdb_taxonomy.tsv",
    "radii/gtdb_radii.tsv",
    "msa/gtdb_r232_bac120.faa",
    "msa/gtdb_r232_ar53.faa",
    "masks/gtdb_r232_bac120.mask",
    "masks/gtdb_r232_ar53.mask",
    "markers/tigrfam/tigrfam.hmm",
    "mrca_red/gtdbtk_r232_bac120.tsv",
    "mrca_red/gtdbtk_r232_ar53.tsv",
)
GTDBTK_REQUIRED_DIRECTORIES = (
    "markers/pfam",
    "pplacer",
    "skani",
    "split",
)
GTDBTK_FATAL_DOWNLOAD_MARKERS = (
    "unauthorized",
    "forbidden",
    "authentication failed",
    "authorization failed",
    "certificate verify failed",
    "certificate verification failed",
    "ssl certificate",
    "tls certificate",
    "permission denied",
    "read-only file system",
    "no space left on device",
    "disk quota exceeded",
    "is a directory",
)
GTDBTK_TRANSIENT_HTTP_STATUS = re.compile(
    r"(?:http[^\n]*\s|error\s+)(429|5\d\d)\b", re.I
)
GTDBTK_FATAL_HTTP_STATUS = re.compile(r"(?:http[^\n]*\s|error\s+)(4(?!29)\d\d)\b", re.I)


def _retryable_gtdbtk_download(error, attempt_output):
    if not isinstance(error, subprocess.CalledProcessError):
        return False
    lowered = attempt_output.lower()
    if any(marker in lowered for marker in GTDBTK_FATAL_DOWNLOAD_MARKERS) or bool(
        GTDBTK_FATAL_HTTP_STATUS.search(attempt_output)
    ):
        return False
    if error.returncode == 4:
        return True
    if error.returncode in {7, 8}:
        return bool(GTDBTK_TRANSIENT_HTTP_STATUS.search(attempt_output))
    return False


def _metadata_release(metadata_path):
    release = None
    try:
        for line in metadata_path.read_text(encoding="utf-8").splitlines():
            key, separator, value = line.partition("=")
            if separator and key.strip() == "VERSION_DATA":
                release = value.strip()
                break
    except (OSError, UnicodeError) as error:
        raise CacheIntegrityError(
            f"Cannot read GTDB-Tk metadata {metadata_path}: {error}"
        ) from error
    return release


def validate_gtdbtk_database(data_dir):
    """Check R232 identity and representative native reference material."""
    data_dir = resolved_path(data_dir)
    if not data_dir.is_dir():
        raise CacheIntegrityError(f"GTDB-Tk database is not a directory: {data_dir}")
    metadata = data_dir / "metadata" / "metadata.txt"
    if not readable_nonempty_file(metadata):
        raise CacheIntegrityError(
            f"GTDB-Tk cache at {data_dir} lacks readable non-empty metadata/metadata.txt."
        )
    release = _metadata_release(metadata)
    if release != GTDBTK_RELEASE:
        raise CacheIntegrityError(
            f"GTDB-Tk 2.7.1 requires VERSION_DATA={GTDBTK_RELEASE}, but selected "
            f"metadata reports {release!r} at {metadata}."
        )
    missing = [
        relative
        for relative in GTDBTK_REQUIRED_FILES
        if not readable_nonempty_file(data_dir / relative)
    ]
    empty_directories = []
    for relative in GTDBTK_REQUIRED_DIRECTORIES:
        directory = data_dir / relative
        try:
            usable = directory.is_dir() and any(directory.iterdir())
        except OSError:
            usable = False
        if not usable:
            empty_directories.append(relative)
    if missing or empty_directories:
        details = []
        if missing:
            details.append("files=" + ",".join(missing))
        if empty_directories:
            details.append("directories=" + ",".join(empty_directories))
        raise CacheIntegrityError(
            f"GTDB-Tk R232 cache at {data_dir} is incomplete ({'; '.join(details)})."
        )
    return {"release": release}


def ensure_gtdbtk_database(data_dir, log_path):
    """Reuse R232 or download/extract/validate it once in isolated staging."""
    data_dir = resolved_path(data_dir)

    def prepare(staging):
        archive = staging / "gtdbtk_r232_data.tar.gz"
        extracted = staging / "database"
        extracted.mkdir()
        retry_download(
            lambda: run_logged(
                ["wget", "--tries=1", GTDBTK_DATABASE_URL, "-O", archive], log_path
            ),
            label="GTDB-Tk R232",
            is_retryable=_retryable_gtdbtk_download,
            append_log=append_log,
            log_path=log_path,
        )
        run_logged(
            [
                "tar",
                "-xzf",
                archive,
                "-C",
                extracted,
                "--strip-components=1",
            ],
            log_path,
        )
        return extracted

    return ensure_fresh_directory_cache(
        data_dir,
        "gtdbtk-r232",
        validate_gtdbtk_database,
        prepare,
        log_path,
        "GTDB-Tk R232",
        (GTDBTK_RELEASE,),
        waiting_reused_only=True,
    )
