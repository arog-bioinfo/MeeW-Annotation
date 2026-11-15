"""MetaEuk/MMseqs2 catalog cache readiness and managed publication."""

from __future__ import annotations

import json
from pathlib import Path

from .cache_runtime import (
    CacheIntegrityError,
    CacheResult,
    append_log,
    atomic_write_json,
    cache_lock,
    publish_directory,
    readable_nonempty_file,
    resolved_path,
    run_logged,
    staging_directory,
)

METAEUK_IDENTITY_NAME = "catalog-identity.json"
METAEUK_IDENTITY_SCHEMA = 1
EXPECTED_MMSEQS_VERSION = "17.b804f"


def metaeuk_native_family(prefix):
    prefix = resolved_path(prefix)
    return (prefix, Path(f"{prefix}.index"), Path(f"{prefix}.dbtype"))


def validate_metaeuk_native_database(prefix):
    """Require the core MMseqs2 data/index/dbtype family, not a lone marker."""
    prefix = resolved_path(prefix)
    missing = [
        path
        for path in metaeuk_native_family(prefix)
        if not readable_nonempty_file(path)
    ]
    if missing:
        raise CacheIntegrityError(
            f"MetaEuk MMseqs2 prefix {prefix} is incomplete; missing, empty, or "
            "unreadable native family members: "
            + ", ".join(str(path) for path in missing)
        )
    return prefix


def identity_path(prefix):
    return resolved_path(prefix).parent / METAEUK_IDENTITY_NAME


def _load_identity(prefix):
    record_path = identity_path(prefix)
    if not record_path.exists():
        return None
    if not readable_nonempty_file(record_path):
        raise CacheIntegrityError(
            f"MetaEuk catalog identity record is empty or unreadable: {record_path}"
        )
    try:
        record = json.loads(record_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise CacheIntegrityError(
            f"MetaEuk catalog identity record is invalid at {record_path}: {error}"
        ) from error
    if not isinstance(record, dict):
        raise CacheIntegrityError(
            f"MetaEuk catalog identity record must be an object: {record_path}"
        )
    return record


def _validate_identity(prefix, catalog, *, required):
    record = _load_identity(prefix)
    if record is None:
        if required:
            raise CacheIntegrityError(
                f"Managed MetaEuk cache at {resolved_path(prefix).parent} has no "
                f"{METAEUK_IDENTITY_NAME}; refusing to guess catalog identity."
            )
        return None
    expected = {
        "schema": METAEUK_IDENTITY_SCHEMA,
        "catalog": catalog,
        "database_prefix": Path(prefix).name,
        "mmseqs_version": EXPECTED_MMSEQS_VERSION,
    }
    mismatches = {
        key: (record.get(key), value)
        for key, value in expected.items()
        if record.get(key) != value
    }
    if mismatches:
        rendered = ", ".join(
            f"{key}={found!r} (expected {expected_value!r})"
            for key, (found, expected_value) in mismatches.items()
        )
        raise CacheIntegrityError(
            f"MetaEuk catalog identity mismatch at {identity_path(prefix)}: {rendered}."
        )
    return record


def inspect_metaeuk_database(prefix, catalog, *, managed):
    prefix = validate_metaeuk_native_database(prefix)
    # A record, wherever present, is authoritative and a known mismatch is fatal.
    record = _validate_identity(prefix, catalog, required=managed)
    return {
        "catalog": catalog,
        "identity": "verified" if record is not None else "asserted_unverified",
    }


def _ready(prefix, catalog, managed):
    try:
        return inspect_metaeuk_database(prefix, catalog, managed=managed)
    except CacheIntegrityError:
        return None


def ensure_metaeuk_database(prefix, catalog, *, managed, threads, log_path):
    """Reuse a complete prefix or prepare only a missing managed catalog leaf."""
    prefix = resolved_path(prefix)
    details = _ready(prefix, catalog, managed)
    if details is not None:
        identity = details["identity"]
        append_log(
            log_path,
            f"MetaEuk catalog {catalog!r} reused read-only from {prefix}; "
            f"identity={identity}.",
        )
        return CacheResult("reused", str(prefix), identity=identity)

    if not managed:
        # Explicit and legacy prefixes are user-owned. Never create or relabel them.
        if any(
            path.exists() or path.is_symlink() for path in metaeuk_native_family(prefix)
        ):
            inspect_metaeuk_database(prefix, catalog, managed=False)
        raise FileNotFoundError(
            f"Configured/legacy MetaEuk prefix {prefix} for catalog {catalog!r} is "
            "missing or incomplete. MeeW will not prepare into an external or legacy "
            "shared prefix; leave metaeuk.db blank and use the managed catalogs cache, "
            "or provision this prefix explicitly."
        )

    leaf = prefix.parent
    if prefix.name != "db":
        raise ValueError(f"Managed MetaEuk prefixes must end in '/db', found {prefix}.")
    with cache_lock(leaf, "metaeuk-catalog"):
        details = _ready(prefix, catalog, managed=True)
        if details is not None:
            append_log(
                log_path,
                f"MetaEuk catalog {catalog!r} became ready while waiting; reused {prefix}.",
            )
            return CacheResult("reused", str(prefix), identity="verified")
        if leaf.exists() or leaf.is_symlink():
            inspect_metaeuk_database(prefix, catalog, managed=True)

        with staging_directory(leaf, "metaeuk-catalog") as staging_container:
            staged_leaf = staging_container / "catalog"
            staged_leaf.mkdir()
            staged_prefix = staged_leaf / "db"
            temporary = staged_leaf / "tmp"
            temporary.mkdir()
            run_logged(
                [
                    "mmseqs",
                    "databases",
                    catalog,
                    staged_prefix,
                    temporary,
                    "--threads",
                    str(threads),
                ],
                log_path,
            )
            validate_metaeuk_native_database(staged_prefix)
            atomic_write_json(
                staged_leaf / METAEUK_IDENTITY_NAME,
                {
                    "schema": METAEUK_IDENTITY_SCHEMA,
                    "catalog": catalog,
                    "database_prefix": "db",
                    "mmseqs_version": EXPECTED_MMSEQS_VERSION,
                    "managed_by": "MeeW-Annotation",
                },
            )
            inspect_metaeuk_database(staged_prefix, catalog, managed=True)
            publish_directory(staged_leaf, leaf)

    inspect_metaeuk_database(prefix, catalog, managed=True)
    append_log(log_path, f"MetaEuk catalog {catalog!r} prepared at {prefix}.")
    return CacheResult("prepared", str(prefix), prepared=(catalog,))
