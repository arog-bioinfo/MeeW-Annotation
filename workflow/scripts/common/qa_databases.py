"""Runtime readiness helpers for Annotation-owned CheckM2 and EukCC adapters."""

from __future__ import annotations

import errno
import http.client
import json
import os
import re
import socket
import ssl
import subprocess
import tarfile
import urllib.error
import urllib.request
from pathlib import Path

from .cache_runtime import (
    CacheIntegrityError,
    CacheResult,
    append_log,
    ensure_fresh_directory_cache,
    readable_nonempty_file,
    resolved_path,
    run_logged,
)
from .download_retries import retry_download

CHECKM2_RELATIVE_DATABASE = Path("CheckM2_database/uniref100.KO.1.dmnd")
EUKCC_DB_URL = (
    "http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/" "eukcc2_db_ver_1.2.tar.gz"
)
EUKCC_ARCHIVE_ROOT = "eukcc2_db_ver_1.2"
CHECKM2_FATAL_DOWNLOAD_MARKERS = (
    "checksum",
    "zenodobackpackmalformedexception",
    "zenodobackpackversionexception",
    "incorrect zenodo backpack version",
    "could not verify successful installation",
    "permission to write",
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
    "404 not found",
    "doi could not be resolved",
)
CHECKM2_TRANSIENT_DOWNLOAD_MARKERS = (
    "connection error:",
    "error during metadata retrieval:",
    "connection timed out",
    "read timed out",
    "temporary failure in name resolution",
    "connection refused",
    "connection reset",
    "connection broken:",
    "incompleteread(",
    "chunkedencodingerror",
    "network is unreachable",
    "no route to host",
    "remotedisconnected",
    "newconnectionerror",
    "failed to establish a new connection",
)
CHECKM2_TRANSIENT_HTTP_STATUS = re.compile(
    r"(?:http[^\n]*\s|status(?: code)?[=: ]+)(429|5\d\d)\b"
    r"|\b(429|5\d\d)\s+(?:server error|internal server error|bad gateway|"
    r"service unavailable|gateway timeout)",
    re.I,
)
CHECKM2_FATAL_HTTP_STATUS = re.compile(
    r"(?:http[^\n]*\s|status(?: code)?[=: ]+)(4(?!29)\d\d)\b"
    r"|\b(4(?!29)\d\d)\s+(?:client error|bad request|unauthorized|forbidden|not found)",
    re.I,
)
TRANSIENT_NETWORK_ERRNOS = {
    errno.ECONNABORTED,
    errno.ECONNREFUSED,
    errno.ECONNRESET,
    errno.EHOSTDOWN,
    errno.EHOSTUNREACH,
    errno.ENETDOWN,
    errno.ENETRESET,
    errno.ENETUNREACH,
    errno.ETIMEDOUT,
}


def _retryable_checkm2_download(error, attempt_output):
    if not isinstance(error, subprocess.CalledProcessError):
        return False
    lowered = attempt_output.lower()
    if any(marker in lowered for marker in CHECKM2_FATAL_DOWNLOAD_MARKERS) or bool(
        CHECKM2_FATAL_HTTP_STATUS.search(attempt_output)
    ):
        return False
    return any(
        marker in lowered for marker in CHECKM2_TRANSIENT_DOWNLOAD_MARKERS
    ) or bool(CHECKM2_TRANSIENT_HTTP_STATUS.search(attempt_output))


def _temporary_url_reason(reason):
    if isinstance(reason, ssl.SSLError):
        return False
    if isinstance(reason, (TimeoutError, socket.timeout)):
        return True
    if isinstance(reason, socket.gaierror):
        return reason.errno == socket.EAI_AGAIN
    if isinstance(reason, OSError):
        return reason.errno in TRANSIENT_NETWORK_ERRNOS
    lowered = str(reason).lower()
    return any(
        marker in lowered
        for marker in (
            "timed out",
            "temporary failure in name resolution",
            "connection refused",
            "connection reset",
            "network is unreachable",
            "no route to host",
        )
    )


def _retryable_eukcc_download(error, _attempt_output):
    if isinstance(error, urllib.error.HTTPError):
        return error.code == 429 or 500 <= error.code <= 599
    if isinstance(error, urllib.error.ContentTooShortError):
        return True
    if isinstance(error, urllib.error.URLError):
        return _temporary_url_reason(error.reason)
    if isinstance(error, ssl.SSLError):
        return False
    if isinstance(
        error,
        (
            TimeoutError,
            ConnectionError,
            http.client.IncompleteRead,
            http.client.RemoteDisconnected,
        ),
    ):
        return True
    if isinstance(error, socket.gaierror):
        return error.errno == socket.EAI_AGAIN
    if isinstance(error, OSError):
        return error.errno in TRANSIENT_NETWORK_ERRNOS
    return False


def selected_checkm2_database(database_path):
    """Return the actual selected .dmnd file for direct-file or root selections."""
    configured = resolved_path(database_path)
    if configured.suffix == ".dmnd" or configured.is_file():
        return configured, None
    return configured / CHECKM2_RELATIVE_DATABASE, configured


def ensure_checkm2_database(database_path, download_database, log_path):
    database_file, root = selected_checkm2_database(database_path)
    if readable_nonempty_file(database_file):
        append_log(log_path, f"CheckM2 database reused read-only from {database_file}.")
        return database_file, CacheResult("reused", str(database_file))
    if root is None:
        raise FileNotFoundError(
            f"Configured direct CheckM2 database file is missing, empty, or unreadable: "
            f"{database_file}. Native download expects a directory root, so MeeW will "
            "not create a directory at this filename."
        )
    if not download_database:
        raise FileNotFoundError(
            f"CheckM2 database file is missing, empty, or unreadable: {database_file}. "
            "Set qa.checkm2.download_database: true to prepare a wholly absent root."
        )

    def validate(selected_root):
        selected_file, _ = selected_checkm2_database(selected_root)
        if not readable_nonempty_file(selected_file):
            if resolved_path(selected_root) != root:
                raise CacheIntegrityError(
                    "CheckM2 download completed without a readable non-empty database "
                    f"at {selected_file}."
                )
            raise CacheIntegrityError(
                f"CheckM2 cache root exists but its selected native database is missing "
                f"or unusable: {database_file}. Refusing automatic repair of {root}."
            )
        return selected_file

    def prepare(staging_container):
        staged_root = staging_container / "database"
        staged_root.mkdir()
        command = [
            "checkm2",
            "database",
            "--download",
            "--path",
            staged_root,
            "--no_write_json_db",
        ]
        retry_download(
            lambda: run_logged(command, log_path),
            label="CheckM2",
            is_retryable=_retryable_checkm2_download,
            append_log=append_log,
            log_path=log_path,
        )
        return staged_root

    result = ensure_fresh_directory_cache(
        root,
        "checkm2",
        validate,
        prepare,
        log_path,
        "CheckM2",
        reported_path=database_file,
    )
    return database_file, CacheResult(result.status, str(database_file))


def _usable_eukcc_path(base, raw):
    path = Path(os.path.expanduser(os.path.expandvars(raw)))
    if not path.is_absolute():
        path = base / path
    if readable_nonempty_file(path):
        return True
    pressed = Path(f"{path}.h3f")
    if readable_nonempty_file(pressed):
        return True
    if path.is_dir():
        try:
            return any(path.iterdir())
        except OSError:
            return False
    return False


def validate_eukcc_database(db_path):
    """Follow EukCC's db_info references and reject empty placeholder directories."""
    db_path = resolved_path(db_path)
    if not db_path.is_dir():
        raise CacheIntegrityError(f"EukCC database is not a directory: {db_path}")
    info_path = db_path / "db_base" / "db_info.json"
    if not readable_nonempty_file(info_path):
        raise CacheIntegrityError(
            f"EukCC cache at {db_path} lacks readable non-empty db_base/db_info.json."
        )
    try:
        info = json.loads(info_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise CacheIntegrityError(
            f"EukCC metadata is invalid at {info_path}: {error}"
        ) from error
    files = info.get("files") if isinstance(info, dict) else None
    if not isinstance(files, dict) or not files:
        raise CacheIntegrityError(
            f"EukCC metadata has no native files map: {info_path}"
        )

    missing = []

    def inspect(mapping, key_path=()):
        for key, value in mapping.items():
            nested_key = (*key_path, str(key))
            if isinstance(value, dict):
                inspect(value, nested_key)
            elif isinstance(value, str):
                if not _usable_eukcc_path(info_path.parent, value):
                    missing.append(".".join(nested_key) + "=" + value)
            else:
                missing.append(".".join(nested_key) + "=<invalid>")

    inspect(files)
    if missing:
        raise CacheIntegrityError(
            f"EukCC cache at {db_path} has missing/empty native material: "
            + ", ".join(missing[:8])
            + (" ..." if len(missing) > 8 else "")
        )
    return info


def _safe_extract(archive_path, destination):
    destination = resolved_path(destination)
    with tarfile.open(archive_path, "r:*") as archive:
        for member in archive.getmembers():
            member_path = destination / member.name
            try:
                member_path.resolve().relative_to(destination)
            except ValueError as error:
                raise CacheIntegrityError(
                    f"Unsafe path in EukCC database archive: {member.name}"
                ) from error
        archive.extractall(destination, filter="data")


def ensure_eukcc_database(db_path, download_database, log_path):
    db_path = resolved_path(db_path)
    if not download_database:
        if db_path.exists() or db_path.is_symlink():
            validate_eukcc_database(db_path)
            append_log(log_path, f"EukCC database reused read-only from {db_path}.")
            return CacheResult("reused", str(db_path))
        raise FileNotFoundError(
            f"EukCC database does not exist at {db_path}; enable "
            "qa.eukcc.download_database to prepare it."
        )

    def prepare(staging):
        archive = staging / Path(EUKCC_DB_URL).name
        append_log(log_path, f"Downloading EukCC database from {EUKCC_DB_URL}.")
        retry_download(
            lambda: urllib.request.urlretrieve(EUKCC_DB_URL, archive),
            label="EukCC",
            is_retryable=_retryable_eukcc_download,
            append_log=append_log,
            log_path=log_path,
        )
        _safe_extract(archive, staging)
        return staging / EUKCC_ARCHIVE_ROOT

    return ensure_fresh_directory_cache(
        db_path,
        "eukcc",
        validate_eukcc_database,
        prepare,
        log_path,
        "EukCC",
    )
