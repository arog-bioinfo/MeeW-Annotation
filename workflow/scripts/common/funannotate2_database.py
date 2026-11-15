"""Funannotate2 26.6.21 selected-database cache lifecycle."""

from __future__ import annotations

import json
import os
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

FUNANNOTATE2_ALL_DATABASES = (
    "merops",
    "uniprot",
    "dbCAN",
    "pfam",
    "go",
    "mibig",
    "interpro",
    "gene2product",
    "mito",
)
FUNANNOTATE2_RECORD_KEYS = ("db", "type", "version", "date", "n_records", "md5")
FUNANNOTATE2_INFO_NAME = "funannotate-db-info.json"


def expand_funannotate2_databases(databases):
    """Match upstream's mixed-``all`` behavior while preserving unique order."""
    requested = tuple(str(database) for database in databases)
    if "all" in requested:
        return FUNANNOTATE2_ALL_DATABASES
    return tuple(dict.fromkeys(requested))


def _load_info(db_dir, *, required):
    info_path = Path(db_dir) / FUNANNOTATE2_INFO_NAME
    if not info_path.exists():
        if required:
            raise CacheIntegrityError(
                f"Funannotate2 cache at {db_dir} is missing {FUNANNOTATE2_INFO_NAME}."
            )
        return {}
    if not readable_nonempty_file(info_path):
        raise CacheIntegrityError(
            f"Funannotate2 metadata is not a readable non-empty file: {info_path}"
        )
    try:
        info = json.loads(info_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise CacheIntegrityError(
            f"Funannotate2 metadata is invalid at {info_path}: {error}"
        ) from error
    if not isinstance(info, dict):
        raise CacheIntegrityError(
            f"Funannotate2 metadata must be a JSON object: {info_path}"
        )
    return info


def _record_database_path(record, db_dir):
    raw = record.get("db")
    if not isinstance(raw, str) or not raw.strip():
        return None
    path = Path(os.path.expanduser(os.path.expandvars(raw)))
    if not path.is_absolute():
        path = Path(db_dir) / path
    return path.resolve()


def validate_funannotate2_entry(db_dir, database, record):
    """Validate one native metadata record and its referenced primary material."""
    if not isinstance(record, dict):
        raise CacheIntegrityError(
            f"Funannotate2 entry {database!r} is not an object in {db_dir}."
        )
    missing = [key for key in FUNANNOTATE2_RECORD_KEYS if key not in record]
    if missing:
        raise CacheIntegrityError(
            f"Funannotate2 entry {database!r} is partial; missing metadata keys: "
            f"{', '.join(missing)}."
        )
    for key in ("type", "version", "date", "md5"):
        if not str(record[key]).strip():
            raise CacheIntegrityError(
                f"Funannotate2 entry {database!r} has blank {key!r} metadata."
            )
    if record["n_records"] in (None, ""):
        raise CacheIntegrityError(
            f"Funannotate2 entry {database!r} has no native record count."
        )
    material = _record_database_path(record, db_dir)
    if material is None or not readable_nonempty_file(material):
        raise CacheIntegrityError(
            f"Funannotate2 entry {database!r} references missing, empty, or unreadable "
            f"material: {record.get('db')!r}."
        )
    return material


def inspect_funannotate2_database(db_dir, databases):
    """Return native metadata and absent requested entries; reject partial entries."""
    db_dir = resolved_path(db_dir)
    required = expand_funannotate2_databases(databases)
    if not db_dir.is_dir():
        raise CacheIntegrityError(
            f"Funannotate2 database cache is not a directory: {db_dir}"
        )
    info = _load_info(db_dir, required=True)
    missing = []
    for database in required:
        if database not in info:
            missing.append(database)
        else:
            validate_funannotate2_entry(db_dir, database, info[database])
    return info, tuple(missing)


def _ready(db_dir, databases):
    try:
        info, missing = inspect_funannotate2_database(db_dir, databases)
    except CacheIntegrityError:
        return None
    return info if not missing else None


def _rebase_staged_records(info, staging_root, final_root):
    rebased = {}
    staging_root = resolved_path(staging_root)
    final_root = resolved_path(final_root)
    for database, record in info.items():
        copied = dict(record)
        material = _record_database_path(copied, staging_root)
        if material is None:
            raise CacheIntegrityError(
                f"Funannotate2 installer emitted an invalid {database!r} db path."
            )
        try:
            relative = material.relative_to(staging_root)
        except ValueError as error:
            raise CacheIntegrityError(
                f"Funannotate2 installer emitted material outside staging for "
                f"{database!r}: {material}"
            ) from error
        copied["db"] = str(final_root / relative)
        rebased[database] = copied
    return rebased


def _install_into_staging(staging, missing, final_root, log_path):
    environment = os.environ.copy()
    environment["FUNANNOTATE2_DB"] = str(staging)
    run_logged(
        ["funannotate2", "install", "-d", *missing],
        log_path,
        env=environment,
    )
    staged_info = _load_info(staging, required=True)
    for database in missing:
        if database not in staged_info:
            raise CacheIntegrityError(
                f"Funannotate2 install completed without requested entry {database!r}."
            )
        validate_funannotate2_entry(staging, database, staged_info[database])
    return _rebase_staged_records(staged_info, staging, final_root)


def _publish_additive(staging, db_dir, current_info, staged_info):
    """Merge isolated new material without replacing existing cache entries."""
    staging = Path(staging)
    db_dir = resolved_path(db_dir)
    additions = [
        child
        for child in staging.iterdir()
        if child.name != FUNANNOTATE2_INFO_NAME
        and not (
            child.name == "pretrained" and child.is_dir() and not any(child.iterdir())
        )
    ]
    collisions = [
        child.name
        for child in additions
        if (db_dir / child.name).exists() or (db_dir / child.name).is_symlink()
    ]
    if collisions:
        raise CacheIntegrityError(
            "Funannotate2 cannot safely add the selected databases because staged "
            "native material collides with existing cache paths: "
            + ", ".join(sorted(collisions))
        )

    moved = []
    try:
        for child in additions:
            destination = db_dir / child.name
            os.replace(child, destination)
            moved.append((destination, child))
        merged = dict(current_info)
        for database, record in staged_info.items():
            if database in merged:
                raise CacheIntegrityError(
                    f"Refusing to replace existing Funannotate2 entry {database!r}."
                )
            merged[database] = record
        atomic_write_json(db_dir / FUNANNOTATE2_INFO_NAME, merged)
    except BaseException:
        for destination, original in reversed(moved):
            if destination.exists() and not original.exists():
                os.replace(destination, original)
        raise


def ensure_funannotate2_database(db_dir, databases, log_path):
    """Reuse a ready selection or safely install only wholly absent entries."""
    db_dir = resolved_path(db_dir)
    required = expand_funannotate2_databases(databases)
    if not required:
        raise ValueError("At least one Funannotate2 database must be selected.")
    if _ready(db_dir, required) is not None:
        append_log(
            log_path,
            f"Funannotate2 database selection reused read-only from {db_dir}: "
            + ", ".join(required),
        )
        return CacheResult("reused", str(db_dir))

    with cache_lock(db_dir, "funannotate2"):
        if _ready(db_dir, required) is not None:
            append_log(
                log_path,
                f"Funannotate2 database selection became ready while waiting; "
                f"reused {db_dir}.",
            )
            return CacheResult("reused", str(db_dir))

        fresh = not db_dir.exists() and not db_dir.is_symlink()
        if fresh:
            current_info = {}
            missing = required
        else:
            if not db_dir.is_dir():
                raise CacheIntegrityError(
                    f"Funannotate2 cache path exists but is not a directory: {db_dir}"
                )
            current_info = _load_info(db_dir, required=False)
            missing = []
            for database in required:
                if database in current_info:
                    validate_funannotate2_entry(
                        db_dir, database, current_info[database]
                    )
                else:
                    missing.append(database)
            missing = tuple(missing)
            if not missing:
                return CacheResult("reused", str(db_dir))

        with staging_directory(db_dir, "funannotate2") as staging_container:
            staging = staging_container / "database"
            staging.mkdir()
            staged_info = _install_into_staging(staging, missing, db_dir, log_path)
            if fresh:
                atomic_write_json(
                    staging / FUNANNOTATE2_INFO_NAME,
                    {database: staged_info[database] for database in missing},
                )
                for database in missing:
                    validate_funannotate2_entry(
                        staging,
                        database,
                        {
                            **staged_info[database],
                            "db": str(
                                staging
                                / Path(staged_info[database]["db"]).relative_to(db_dir)
                            ),
                        },
                    )
                publish_directory(staging, db_dir)
            else:
                _publish_additive(
                    staging,
                    db_dir,
                    current_info,
                    {database: staged_info[database] for database in missing},
                )

        inspect_funannotate2_database(db_dir, required)

    append_log(
        log_path,
        f"Funannotate2 prepared missing database entries at {db_dir}: "
        + ", ".join(missing),
    )
    return CacheResult("prepared", str(db_dir), prepared=tuple(missing))
