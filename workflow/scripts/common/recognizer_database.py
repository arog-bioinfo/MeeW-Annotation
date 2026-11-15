"""Readiness and lock assessment for native reCOGnizer 1.12.2 resources."""

from __future__ import annotations

import os
import shlex
import shutil
from dataclasses import dataclass
from pathlib import Path

from .cache_runtime import CacheIntegrityError, readable_nonempty_file, resolved_path
from .native_execution import exclusive_lock, shared_lock

RECOGNIZER_COMPLETION_INDICATOR = "recognizer_dwnl.timestamp"
RECOGNIZER_INDEX_EXTENSIONS = (
    "aux",
    "freq",
    "loo",
    "pdb",
    "phr",
    "pin",
    "pos",
    "pot",
    "psq",
    "ptf",
    "pto",
    "rps",
)
RECOGNIZER_PREFIXES = {
    "NCBI_Curated": ("Cdd_NCBI", ("cd", "sd")),
    "Pfam": ("Pfam", ("pfam",)),
    "SMART": ("Smart", ("smart",)),
    "KOG": ("Kog", ("KOG",)),
    "COG": ("Cog", ("COG",)),
    "PRK": ("Prk", ("CHL", "MTH", "NF", "PHA", "PLN", "PRK", "PTZ")),
    "TIGR": ("Tigr", ("TIGR",)),
}
RECOGNIZER_COMMON_RELATIONAL_FILES = (
    "bitscore_specific.txt",
    "cddannot.dat",
    "cddannot_generic.dat",
    "cddid.tbl",
    "cdtrack.txt",
    "family_superfamily_links",
    "hmm_PGAP.tsv",
)
RECOGNIZER_EMPTY_COMMON_RELATIONAL_FILE = "cdtrack.txt"
RECOGNIZER_DATABASE_RELATIONAL_FILES = {
    "SMART": ("descriptions.pl",),
    "KOG": ("kog",),
    "COG": ("cog-20.def.tab",),
}
RECOGNIZER_COG2EC_SOURCES = (
    "eggnog4.protein_id_conversion.tsv",
    "NOG.members.tsv",
)
RECOGNIZER_COG2KO_SOURCES = (
    "COG.mappings.v11.0.txt",
    "protein.info.v11.0.txt",
)
RECOGNIZER_COG_MAP_SOURCES = (*RECOGNIZER_COG2EC_SOURCES, *RECOGNIZER_COG2KO_SOURCES)
RECOGNIZER_COG_MAP_SOURCE_GROUPS = {
    "cog2ec_recognizer.tsv": RECOGNIZER_COG2EC_SOURCES,
    "cog2ko_recognizer.tsv": RECOGNIZER_COG2KO_SOURCES,
}


@dataclass(frozen=True)
class RecognizerAssessment:
    ready: bool
    may_mutate: bool
    locks: tuple
    reasons: tuple[str, ...]
    resources_root: Path
    native_source: Path


def _database_names(databases):
    if isinstance(databases, str):
        return tuple(name for name in databases.split(",") if name)
    return tuple(str(name) for name in databases)


def _native_source(executable=None):
    selected = executable or shutil.which("recognizer")
    if not selected:
        raise CacheIntegrityError(
            "Cannot locate the native reCOGnizer executable in PATH, so its "
            "installation-level mapping resources cannot be protected."
        )
    source = Path(selected).expanduser().resolve()
    if not source.is_file():
        raise CacheIntegrityError(
            f"Native reCOGnizer entrypoint does not resolve to a file: {source}"
        )
    return source


def recognizer_native_assets(executable=None):
    source = _native_source(executable)
    directory = source.parent
    return (
        source,
        directory / "fun.tsv",
        (
            directory / "cog2ec_recognizer.tsv",
            directory / "cog2ko_recognizer.tsv",
        ),
    )


def _readable_file(path, allow_empty=False):
    if not allow_empty:
        return readable_nonempty_file(path)
    path = Path(path)
    try:
        return path.is_file() and os.access(path, os.R_OK)
    except OSError:
        return False


def _require_file(path, description, allow_empty=False):
    if not _readable_file(path, allow_empty):
        failure = (
            "missing or unreadable" if allow_empty else "missing, empty, or unreadable"
        )
        raise CacheIntegrityError(
            f"Incomplete reCOGnizer cache: required {description} is {failure} "
            f"at {path}. The guard will not claim this "
            "selection is ready or repair it automatically."
        )


def _file_ready(path, description, definitive, allow_empty=False):
    ready = _readable_file(path, allow_empty)
    if not ready and definitive:
        _require_file(path, description, allow_empty)
    return ready


def _has_payload(root, definitive):
    if not root.exists():
        return False
    try:
        return any(path.is_file() or path.is_symlink() for path in root.rglob("*"))
    except OSError as error:
        if definitive:
            raise CacheIntegrityError(
                f"Cannot inspect reCOGnizer resource directory {root}: {error}"
            ) from error
        return True


def _extra_may_write(extra):
    tokens = shlex.split("" if extra is None else str(extra), posix=True)
    harmless_flags = {
        "--debug",
        "--keep-spaces",
        "--no-blast-info",
        "--no-output-sequences",
        "--output-rpsbproc-cols",
        "--quiet",
    }
    harmless_values = {
        "--evalue",
        "--max-target-seqs",
    }
    short_values = {"-mts"}
    index = 0
    while index < len(tokens):
        token = tokens[index]
        option = token.split("=", 1)[0]
        if option in harmless_flags:
            index += 1
            continue
        if option in harmless_values or option in short_values:
            if "=" not in token:
                index += 1
            index += 1
            continue
        # Taxonomic/test modes and shell constructs are intentionally
        # conservative: native 1.12.2 can create additional shared indices.
        return bool(tokens)
    return False


def recognizer_material_candidates(resources_root, databases):
    """Return selected native material paths without reading or creating them."""
    root = resolved_path(resources_root or "~/recognizer_resources")
    names = _database_names(databases)
    candidates = [root / RECOGNIZER_COMPLETION_INDICATOR]
    candidates.extend(root / name for name in RECOGNIZER_COMMON_RELATIONAL_FILES)
    for database in names:
        candidates.extend(
            root / name
            for name in RECOGNIZER_DATABASE_RELATIONAL_FILES.get(database, ())
        )
        prefix = RECOGNIZER_PREFIXES.get(database)
        if prefix:
            candidates.extend(
                root / "dbs" / f"{prefix[0]}.{extension}"
                for extension in RECOGNIZER_INDEX_EXTENSIONS
            )
    return tuple(candidates)


def assess_recognizer_database(
    resources_root,
    databases,
    extra="",
    *,
    executable=None,
    definitive=True,
):
    """Classify one original native invocation as verified-reader or writer."""
    root = resolved_path(resources_root or "~/recognizer_resources")
    names = _database_names(databases)
    source, fun_table, cog_maps = recognizer_native_assets(executable)
    reasons = []
    native_writes = not _file_ready(fun_table, "installation-level fun.tsv", definitive)
    if native_writes:
        reasons.append("unusable installation-level fun.tsv")

    root_is_directory = not root.exists() or root.is_dir()
    if not root_is_directory:
        if definitive:
            raise CacheIntegrityError(
                f"reCOGnizer resources path exists but is not a directory: {root}"
            )
        reasons.append("unusable resources path")

    timestamp = root / RECOGNIZER_COMPLETION_INDICATOR
    timestamp_ready = root_is_directory and readable_nonempty_file(timestamp)
    resource_writes = not timestamp_ready or not root_is_directory
    if root_is_directory and timestamp.exists() and not timestamp_ready:
        if definitive:
            _require_file(timestamp, "native completion indicator")
        reasons.append("unusable native completion indicator")
    if root_is_directory and not timestamp_ready and _has_payload(root, definitive):
        if definitive:
            raise CacheIntegrityError(
                f"reCOGnizer resources at {root} contain material but lack a usable "
                f"{RECOGNIZER_COMPLETION_INDICATOR}. This is an incomplete native "
                "initialization; preserve it for diagnosis and explicitly recover it "
                "before retrying."
            )
        reasons.append("unverified in-progress or incomplete initialization")

    if timestamp_ready:
        for name in RECOGNIZER_COMMON_RELATIONAL_FILES:
            if not _file_ready(
                root / name,
                f"relational resource {name}",
                definitive,
                allow_empty=name == RECOGNIZER_EMPTY_COMMON_RELATIONAL_FILE,
            ):
                resource_writes = True
                reasons.append(f"unusable relational resource {name}")
        for database in names:
            for name in RECOGNIZER_DATABASE_RELATIONAL_FILES.get(database, ()):
                if not _file_ready(
                    root / name,
                    f"{database} relational resource {name}",
                    definitive,
                ):
                    resource_writes = True
                    reasons.append(f"unusable {database} relational resource {name}")
            selected = RECOGNIZER_PREFIXES.get(database)
            if selected is None:
                resource_writes = True
                reasons.append(f"unproven database selection {database}")
                continue
            prefix, smp_prefixes = selected
            smps = [
                path
                for smp_prefix in smp_prefixes
                for path in (root / "smps").glob(f"{smp_prefix}*.smp")
            ]
            if not smps or not all(readable_nonempty_file(path) for path in smps):
                if definitive:
                    raise CacheIntegrityError(
                        f"Incomplete reCOGnizer cache: selected {database} has no "
                        f"usable SMP source family in {root / 'smps'}."
                    )
                resource_writes = True
                reasons.append(f"unusable {database} SMP source family")
            family = tuple(
                root / "dbs" / f"{prefix}.{extension}"
                for extension in RECOGNIZER_INDEX_EXTENSIONS
            )
            existing = tuple(path for path in family if path.exists())
            if existing and not all(readable_nonempty_file(path) for path in family):
                missing = [
                    str(path) for path in family if not readable_nonempty_file(path)
                ]
                if definitive:
                    raise CacheIntegrityError(
                        f"Incomplete reCOGnizer {database} RPS index family for {prefix}; "
                        f"missing, empty, or unreadable members: {', '.join(missing)}. "
                        "The guard will not expose this partial family to readers."
                    )
                resource_writes = True
                reasons.append(f"unverified partial {database} RPS index")
            if not existing:
                resource_writes = True
                reasons.append(f"missing lazy {database} RPS index")

    if "COG" in names:
        for path in cog_maps:
            map_ready = readable_nonempty_file(path)
            if path.exists() and not map_ready:
                if definitive:
                    raise CacheIntegrityError(
                        "Incomplete installation-level reCOGnizer COG mapping file: "
                        f"{path}"
                    )
                reasons.append(f"unusable native COG mapping {path.name}")
            if not map_ready:
                native_writes = True
                resource_writes = True
                reasons.append(f"missing native derived mapping {path.name}")
                if timestamp_ready:
                    for name in RECOGNIZER_COG_MAP_SOURCE_GROUPS[path.name]:
                        if not _file_ready(
                            root / name,
                            f"{path.name} mapping source {name}",
                            definitive,
                        ):
                            resource_writes = True
                            reasons.append(f"unusable mapping source {name}")

    if _extra_may_write(extra):
        resource_writes = True
        native_writes = True
        reasons.append("advanced or unproven extra options")

    resource_request = (
        exclusive_lock(root, "recognizer-resources", "; ".join(reasons))
        if resource_writes
        else shared_lock(root, "recognizer-resources", "verified native reader")
    )
    native_request = (
        exclusive_lock(source, "recognizer-native", "; ".join(reasons))
        if native_writes
        else shared_lock(source, "recognizer-native", "installed mappings reader")
    )
    may_mutate = resource_writes or native_writes
    return RecognizerAssessment(
        ready=not may_mutate,
        may_mutate=may_mutate,
        locks=(resource_request, native_request),
        reasons=tuple(reasons),
        resources_root=root,
        native_source=source,
    )
