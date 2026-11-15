"""Bakta 1.12.0 database readiness and first-install lifecycle."""

from __future__ import annotations

import json
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

BAKTA_DB_SCHEMA_VERSION = 6
BAKTA_COMMON_FILES = (
    "antifam.h3f",
    "antifam.h3i",
    "antifam.h3m",
    "antifam.h3p",
    "bakta.db",
    "expert-protein-sequences.dmnd",
    "ncRNA-genes.i1f",
    "ncRNA-genes.i1i",
    "ncRNA-genes.i1m",
    "ncRNA-genes.i1p",
    "ncRNA-regions.i1f",
    "ncRNA-regions.i1i",
    "ncRNA-regions.i1m",
    "ncRNA-regions.i1p",
    "oric.fna",
    "orit.fna",
    "pfam.h3f",
    "pfam.h3i",
    "pfam.h3m",
    "pfam.h3p",
    "rfam-go.tsv",
    "rRNA.i1f",
    "rRNA.i1i",
    "rRNA.i1m",
    "rRNA.i1p",
    "sorf.dmnd",
)
BAKTA_DOWNLOAD_FAILURE = "ERROR: Could not download file from Zenodo!"
BAKTA_FATAL_DOWNLOAD_MARKERS = (
    "checksum",
    "md5 mismatch",
    "could not extract",
    "error extracting",
    "amrfinderplus",
    "permission denied",
    "read-only file system",
    "no space left on device",
    "disk quota exceeded",
    "certificate verify failed",
    "certificate_verify_failed",
    "ssl certificate",
    "tls certificate",
    "unauthorized",
    "forbidden",
    "authentication failed",
    "authorization failed",
)
BAKTA_FATAL_HTTP_STATUS = re.compile(
    r"\b(?:http(?:/\d(?:\.\d)?)?|http error|error)\s+(401|403|404)\b"
    r"|\bstatus(?: code)?\s*[=:]\s*(401|403|404)\b"
    r"|\b(401|403|404)\s+(?:client error|unauthorized|forbidden|not found)\b",
    re.I,
)


def _retryable_bakta_download(error, attempt_output):
    if not isinstance(error, subprocess.CalledProcessError):
        return False
    lowered = attempt_output.lower()
    if any(marker in lowered for marker in BAKTA_FATAL_DOWNLOAD_MARKERS) or bool(
        BAKTA_FATAL_HTTP_STATUS.search(attempt_output)
    ):
        return False
    return BAKTA_DOWNLOAD_FAILURE.lower() in lowered


def validate_bakta_database(db_dir, selected_type):
    """Validate the pinned native metadata and Bakta's required file set."""
    db_dir = resolved_path(db_dir)
    if selected_type not in {"full", "light"}:
        raise ValueError("Bakta database type must be 'full' or 'light'.")
    if not db_dir.is_dir():
        raise CacheIntegrityError(f"Bakta database is not a directory: {db_dir}")

    version_path = db_dir / "version.json"
    if not readable_nonempty_file(version_path):
        raise CacheIntegrityError(
            f"Bakta cache at {db_dir} is missing readable non-empty version.json."
        )
    try:
        metadata = json.loads(version_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise CacheIntegrityError(
            f"Bakta cache has invalid version metadata at {version_path}: {error}"
        ) from error

    missing_keys = [
        key for key in ("date", "major", "minor", "type") if key not in metadata
    ]
    if missing_keys:
        raise CacheIntegrityError(
            f"Bakta version metadata is missing keys: {', '.join(missing_keys)}."
        )
    if metadata["major"] != BAKTA_DB_SCHEMA_VERSION:
        raise CacheIntegrityError(
            f"Bakta 1.12.0 requires database schema {BAKTA_DB_SCHEMA_VERSION}.x, "
            f"but {metadata['major']!r}.{metadata['minor']!r} is selected at {db_dir}."
        )
    if not isinstance(metadata["minor"], int) or not str(metadata["date"]).strip():
        raise CacheIntegrityError(
            f"Bakta version metadata has invalid date/minor values at {version_path}."
        )
    if metadata["type"] != selected_type:
        raise CacheIntegrityError(
            f"Configured Bakta database type is {selected_type!r}, but cache metadata "
            f"reports {metadata['type']!r} at {db_dir}."
        )

    required = (
        *BAKTA_COMMON_FILES,
        "psc.dmnd" if selected_type == "full" else "pscc.dmnd",
    )
    missing = [name for name in required if not readable_nonempty_file(db_dir / name)]
    if missing:
        preview = ", ".join(missing[:5])
        suffix = " ..." if len(missing) > 5 else ""
        raise CacheIntegrityError(
            f"Bakta {selected_type} cache at {db_dir} is incomplete; missing or "
            f"unreadable native material: {preview}{suffix}"
        )
    return metadata


def ensure_bakta_database(db_dir, selected_type, log_path):
    """Reuse a ready cache or install once into a missing selected directory."""
    db_dir = resolved_path(db_dir)

    def prepare(staging):
        command = [
            "bakta_db",
            "download",
            "--type",
            selected_type,
            "--output",
            staging,
        ]
        # Bakta 1.12.0 opens its archive with ``wb``. Every retry therefore
        # restarts that native transfer from byte zero; this is not resumability.
        retry_download(
            lambda: run_logged(command, log_path),
            label=f"Bakta {selected_type}",
            is_retryable=_retryable_bakta_download,
            append_log=append_log,
            log_path=log_path,
        )
        return staging / ("db" if selected_type == "full" else "db-light")

    return ensure_fresh_directory_cache(
        db_dir,
        "bakta",
        lambda path: validate_bakta_database(path, selected_type),
        prepare,
        log_path,
        f"Bakta {selected_type}",
        (selected_type,),
        detail=lambda metadata: f"schema {metadata['major']}.{metadata['minor']}",
    )
