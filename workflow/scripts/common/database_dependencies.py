"""Lightweight read-only database inputs for Snakemake provenance."""

from __future__ import annotations

import json
import os
from pathlib import Path

from .bakta_database import BAKTA_COMMON_FILES
from .funannotate2_database import (
    FUNANNOTATE2_INFO_NAME,
    expand_funannotate2_databases,
)
from .gtdbtk_database import GTDBTK_REQUIRED_FILES
from .metaeuk_database import METAEUK_IDENTITY_NAME, metaeuk_native_family
from .recognizer_database import recognizer_material_candidates
from .upimapi_database import upimapi_material_candidates


def _available_files(paths):
    """Return existing file dependencies without asserting runtime readiness."""
    available = []
    seen = set()
    for candidate in paths:
        path = Path(candidate).expanduser().resolve()
        if path in seen:
            continue
        seen.add(path)
        if path.is_file():
            available.append(str(path))
    return available


def bakta_database_inputs(db_dir, selected_type):
    """Expose available pinned Bakta identity/material files as read-only inputs."""
    db_dir = Path(db_dir).expanduser().resolve()
    type_file = "psc.dmnd" if selected_type == "full" else "pscc.dmnd"
    return _available_files(
        db_dir / name for name in ("version.json", *BAKTA_COMMON_FILES, type_file)
    )


def gtdbtk_database_inputs(data_dir):
    """Expose available representative R232 identity/material files."""
    data_dir = Path(data_dir).expanduser().resolve()
    return _available_files(data_dir / relative for relative in GTDBTK_REQUIRED_FILES)


def funannotate2_database_inputs(db_dir, databases):
    """Expose native metadata and available selected primary database material."""
    db_dir = Path(db_dir).expanduser().resolve()
    info_path = db_dir / FUNANNOTATE2_INFO_NAME
    candidates = [info_path]
    if info_path.is_file():
        try:
            info = json.loads(info_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError):
            info = None
        if isinstance(info, dict):
            for database in expand_funannotate2_databases(databases):
                record = info.get(database)
                raw = record.get("db") if isinstance(record, dict) else None
                if isinstance(raw, str) and raw.strip():
                    material = Path(os.path.expanduser(os.path.expandvars(raw)))
                    if not material.is_absolute():
                        material = db_dir / material
                    candidates.append(material)
    return _available_files(candidates)


def metaeuk_database_inputs(prefix):
    """Expose available native MMseqs2 family and managed identity files."""
    prefix = Path(prefix).expanduser().resolve()
    return _available_files(
        (*metaeuk_native_family(prefix), prefix.parent / METAEUK_IDENTITY_NAME)
    )


def recognizer_database_inputs(resources_dir, databases):
    """Expose available native completion, relational, and selected index files."""
    return _available_files(recognizer_material_candidates(resources_dir, databases))


def upimapi_database_inputs(resources_dir, db, db_custom):
    """Expose the available selected FASTA/index pair without probing or writing."""
    return _available_files(upimapi_material_candidates(resources_dir, db, db_custom))
