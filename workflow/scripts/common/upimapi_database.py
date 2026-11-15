"""Readiness and lock assessment for native UPIMAPI 1.13.3 resources."""

from __future__ import annotations

import os
import shlex
import struct
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from .cache_runtime import CacheIntegrityError, readable_nonempty_file, resolved_path
from .native_execution import exclusive_lock, shared_lock

UPIMAPI_BUILTIN_FASTAS = {
    "uniprot": "uniprot.fasta",
    "swissprot": "uniprot_sprot.fasta",
    "taxids": "taxids_database.fasta",
}
DIAMOND_MAGIC = bytes.fromhex("6d18ee15a4f84a02")
# Serialized sizes from DIAMOND 2.2.6 ReferenceHeader, its dynamic-record
# length prefix, and SequenceFile::SeqInfo. Bounds below are derived from the
# index's own offsets and counts rather than an arbitrary database file size.
DIAMOND_REFERENCE_HEADER_SIZE = 40
DIAMOND_DYNAMIC_HEADER_PREFIX_SIZE = 8
DIAMOND_SEQUENCE_INFO_SIZE = 16


@dataclass(frozen=True)
class UpimapiAssessment:
    ready: bool
    may_mutate: bool
    locks: tuple
    reasons: tuple[str, ...]
    resources_root: Path
    reference: Path
    diamond_index: Path


def _effective_resources_root(resources_dir):
    return resolved_path(resources_dir or "~/upimapi_resources")


def selected_upimapi_reference(resources_dir, db, db_custom):
    """Mirror native path selection without changing the original CLI string."""
    resources = Path(os.path.expanduser(str(resources_dir or "~/upimapi_resources")))
    if db_custom:
        selected = str(db_custom)
    else:
        database = str(db or "uniprot")
        selected = str(resources / UPIMAPI_BUILTIN_FASTAS.get(database, database))
    if selected.endswith(".dmnd"):
        index = selected
    else:
        index = f"{'.'.join(selected.split('.')[:-1])}.dmnd"
    return resolved_path(selected), resolved_path(index)


def _validate_fasta(path, *, definitive=True):
    if not path.exists():
        return False
    if not readable_nonempty_file(path):
        if definitive:
            raise CacheIntegrityError(
                f"Selected UPIMAPI FASTA is empty, unreadable, or not a regular file: {path}"
            )
        return False
    try:
        with path.open("rb") as handle:
            start = handle.read(65536).lstrip()
    except OSError as error:
        if definitive:
            raise CacheIntegrityError(
                f"Cannot inspect UPIMAPI FASTA {path}: {error}"
            ) from error
        return False
    if not start.startswith(b">"):
        if definitive:
            raise CacheIntegrityError(
                f"Selected UPIMAPI FASTA does not begin with a FASTA header: {path}"
            )
        return False
    return True


def _diamond_layout(path):
    """Inspect fixed metadata and bounded first/terminal layout records."""
    size = path.stat().st_size
    prefix_size = DIAMOND_REFERENCE_HEADER_SIZE + DIAMOND_DYNAMIC_HEADER_PREFIX_SIZE
    with path.open("rb") as handle:
        prefix = handle.read(prefix_size)
        if len(prefix) != prefix_size or not prefix.startswith(DIAMOND_MAGIC):
            return None
        build, format_version, sequence_count, letter_count, positions_offset = (
            struct.unpack(
                "<IIQQQ",
                prefix[len(DIAMOND_MAGIC) : DIAMOND_REFERENCE_HEADER_SIZE],
            )
        )
        (dynamic_header_size,) = struct.unpack(
            "<Q", prefix[DIAMOND_REFERENCE_HEADER_SIZE:prefix_size]
        )
        data_offset = prefix_size + dynamic_header_size
        positions_end = (
            positions_offset + (sequence_count + 1) * DIAMOND_SEQUENCE_INFO_SIZE
        )
        if not (
            sequence_count > 0
            and letter_count > 0
            and data_offset <= positions_offset <= size
            and positions_end <= size
        ):
            return None

        handle.seek(positions_offset)
        first_raw = handle.read(DIAMOND_SEQUENCE_INFO_SIZE)
        second_raw = handle.read(DIAMOND_SEQUENCE_INFO_SIZE)
        handle.seek(positions_offset + sequence_count * DIAMOND_SEQUENCE_INFO_SIZE)
        terminal_raw = handle.read(DIAMOND_SEQUENCE_INFO_SIZE)
        if any(
            len(value) != DIAMOND_SEQUENCE_INFO_SIZE
            for value in (first_raw, second_raw, terminal_raw)
        ):
            return None
        first_position, first_length, first_reserved = struct.unpack("<QII", first_raw)
        second_position, second_length, second_reserved = struct.unpack(
            "<QII", second_raw
        )
        terminal_position, terminal_length, terminal_reserved = struct.unpack(
            "<QII", terminal_raw
        )
        if not (
            first_position == data_offset
            and 0 < first_length <= letter_count
            and first_reserved == 0
            and first_position + first_length + 3 <= second_position
            and second_position <= positions_offset
            and second_reserved == 0
            and terminal_position == positions_offset
            and terminal_length == 0
            and terminal_reserved == 0
            and (sequence_count > 1 or first_length == letter_count)
            and (sequence_count == 1 or second_length > 0)
        ):
            return None

        handle.seek(first_position)
        record_start = handle.read(1)
        handle.seek(first_position + first_length + 1)
        sequence_end = handle.read(1)
        handle.seek(second_position - 1)
        title_end = handle.read(1)
        if record_start != b"\xff" or sequence_end != b"\xff" or title_end != b"\x00":
            return None

    return {
        "build": build,
        "format_version": format_version,
        "sequence_count": sequence_count,
        "letter_count": letter_count,
    }


def _probe_diamond(path, probe):
    if probe is None:
        raise CacheIntegrityError(
            "Definitive UPIMAPI DIAMOND validation requires a supervised native "
            "probe while the selected material locks are held."
        )
    with tempfile.TemporaryDirectory(
        prefix="meew-diamond-readonly-", dir="/tmp"
    ) as scratch:
        output = Path(scratch) / "first-record.faa"
        commands = (
            ("metadata", ("diamond", "dbinfo", "--db", str(path))),
            (
                "first record",
                (
                    "diamond",
                    "getseq",
                    "--db",
                    str(path),
                    "--seq",
                    "1",
                    "--out",
                    str(output),
                ),
            ),
        )
        for label, command in commands:
            try:
                probe(command, cwd=scratch)
            except (OSError, subprocess.CalledProcessError) as error:
                raise CacheIntegrityError(
                    f"Pinned DIAMOND could not read {label} from UPIMAPI index "
                    f"{path}: {error}"
                ) from error
        try:
            with output.open("rb") as handle:
                start = handle.read(65536)
        except OSError as error:
            raise CacheIntegrityError(
                f"Pinned DIAMOND did not produce a readable first record for "
                f"UPIMAPI index {path}: {error}"
            ) from error
        if not start.startswith(b">"):
            raise CacheIntegrityError(
                f"Pinned DIAMOND did not extract a FASTA record from UPIMAPI "
                f"index {path}."
            )


def _validate_diamond(path, *, definitive=True, probe=None):
    if not path.exists():
        return False
    if not readable_nonempty_file(path):
        if definitive:
            raise CacheIntegrityError(
                f"Selected UPIMAPI DIAMOND index is empty, unreadable, or not a file: {path}"
            )
        return False
    try:
        layout = _diamond_layout(path)
    except OSError as error:
        if definitive:
            raise CacheIntegrityError(
                f"Cannot inspect UPIMAPI DIAMOND index {path}: {error}"
            ) from error
        return False
    if layout is None:
        if definitive:
            raise CacheIntegrityError(
                f"Selected UPIMAPI DIAMOND index has an incomplete or invalid "
                f"DIAMOND header/body: {path}. The guard will not treat an "
                "existing truncated or corrupt index as complete."
            )
        return False
    if definitive:
        _probe_diamond(path, probe)
    return True


def _extra_may_write(extra):
    tokens = shlex.split("" if extra is None else str(extra), posix=True)
    harmless_flags = {"--blast", "--show-available-fields"}
    harmless_values = {
        "--bitscore",
        "--block-size",
        "--columns",
        "--diamond-mode",
        "--evalue",
        "--full-id",
        "--index-chunks",
        "--max-memory",
        "--max-tries",
        "--mirror",
        "--pident",
        "--sleep",
        "--step",
    }
    short_values = {"-b", "-c", "-cols"}
    mutating = {"--local-id-mapping", "--taxids"}
    index = 0
    while index < len(tokens):
        token = tokens[index]
        option = token.split("=", 1)[0]
        if option in mutating:
            return True
        if option in harmless_flags:
            index += 1
            continue
        if option in harmless_values or option in short_values:
            if "=" not in token:
                index += 1
            index += 1
            continue
        return bool(tokens)
    return False


def upimapi_material_candidates(resources_dir, db, db_custom):
    """Return selected native material paths without reading or creating them."""
    reference, index = selected_upimapi_reference(resources_dir, db, db_custom)
    return (index,) if reference == index else (reference, index)


def assess_upimapi_database(
    resources_dir,
    db,
    db_custom,
    *,
    skip_db_check,
    extra="",
    definitive=True,
    diamond_probe=None,
):
    """Classify the exact original native invocation's shared-resource access."""
    root = _effective_resources_root(resources_dir)
    root_usable = not root.exists() or root.is_dir()
    if not root_usable and definitive:
        raise CacheIntegrityError(
            f"UPIMAPI resources path exists but is not a directory: {root}"
        )
    reference, index = selected_upimapi_reference(resources_dir, db, db_custom)
    custom = bool(db_custom)
    direct_index = reference == index
    fasta_ready = (
        True if direct_index else _validate_fasta(reference, definitive=definitive)
    )
    index_ready = _validate_diamond(
        index,
        definitive=definitive,
        probe=diamond_probe,
    )
    if not fasta_ready and index_ready:
        if definitive:
            raise CacheIntegrityError(
                f"UPIMAPI index {index} exists without its selected FASTA {reference}; "
                "the guard cannot verify this nonempty incomplete pair."
            )

    reasons = []
    root_writes = not root.exists() or not root_usable
    source_writes = not fasta_ready
    index_writes = not index_ready
    if not fasta_ready:
        reasons.append("missing selected FASTA")
    if not index_ready:
        reasons.append("missing DIAMOND index")
    if not skip_db_check:
        # For built-ins this can download/overwrite. Custom handling is kept
        # native and is conservatively exclusive rather than reinterpreted.
        source_writes = True
        index_writes = True
        if not custom:
            root_writes = True
        reasons.append("native database check enabled")
    advanced = _extra_may_write(extra) or (
        not custom and str(db or "uniprot") == "taxids"
    )
    if advanced:
        root_writes = True
        source_writes = True
        index_writes = True
        reasons.append("taxonomic, local-ID, or unproven extra options")

    for material, may_write in (
        (reference, source_writes),
        (index, index_writes),
    ):
        try:
            material.relative_to(root)
            inside_root = True
        except ValueError:
            inside_root = False
        if may_write and inside_root:
            root_writes = True
    if not custom and (source_writes or index_writes):
        root_writes = True

    resource_request = (
        exclusive_lock(root, "upimapi-resources", "; ".join(reasons))
        if root_writes
        else shared_lock(root, "upimapi-resources", "verified native reader")
    )
    locks = [resource_request]
    if not direct_index:
        source_request = (
            exclusive_lock(reference, "upimapi-material", "; ".join(reasons))
            if source_writes
            else shared_lock(reference, "upimapi-material", "verified reference reader")
        )
        locks.append(source_request)
    target_request = (
        exclusive_lock(index, "upimapi-material", "; ".join(reasons))
        if index_writes
        else shared_lock(index, "upimapi-material", "verified index reader")
    )
    locks.append(target_request)

    may_mutate = root_writes or source_writes or index_writes
    return UpimapiAssessment(
        ready=definitive and not may_mutate,
        may_mutate=may_mutate,
        locks=tuple(locks),
        reasons=tuple(reasons),
        resources_root=root,
        reference=reference,
        diamond_index=index,
    )
