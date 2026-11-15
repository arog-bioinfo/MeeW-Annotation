"""Normalize database selections and derive their persistent resource paths."""

from dataclasses import dataclass
from pathlib import Path
from urllib.parse import quote

from .paths import expand_path, resolve_resource_path

DEFAULT_METAEUK_DATABASE = "UniRef90"
DEFAULT_BAKTA_TYPE = "full"
BAKTA_TYPES = {"full", "light"}


@dataclass(frozen=True)
class MetaEukDatabaseLocation:
    catalog: str
    prefix: str
    managed: bool


def normalize_metaeuk_database(value):
    """Return the selected MMseqs2 catalog name, using UniRef90 when blank."""
    if value is None:
        return DEFAULT_METAEUK_DATABASE
    if not isinstance(value, str):
        raise TypeError("config.metaeuk.database must be a string when provided")
    selected = value.strip()
    return selected or DEFAULT_METAEUK_DATABASE


def metaeuk_database_child(selected_database):
    """Encode one selected catalog name as a safe, deterministic cache child."""
    if selected_database == DEFAULT_METAEUK_DATABASE:
        return "uniref90"
    return f"catalog-{quote(selected_database, safe='')}"


def resolve_metaeuk_database_location(database, explicit_db, cache_dir):
    """Resolve explicit, legacy, or collision-safe managed catalog placement."""
    selected = normalize_metaeuk_database(database)
    explicit = expand_path(explicit_db)
    if explicit:
        return MetaEukDatabaseLocation(selected, explicit, False)

    cache_dir = Path(cache_dir)
    child = metaeuk_database_child(selected)
    managed = cache_dir / "catalogs" / child / "db"
    legacy = cache_dir / child
    # Prefer any new managed leaf, including a partial one that must fail safely,
    # then preserve complete historical flat-prefix selection without relabeling it.
    if managed.parent.exists():
        return MetaEukDatabaseLocation(selected, str(managed), True)
    if any(
        path.exists()
        for path in (legacy, Path(f"{legacy}.index"), Path(f"{legacy}.dbtype"))
    ):
        return MetaEukDatabaseLocation(selected, str(legacy), False)
    return MetaEukDatabaseLocation(selected, str(managed), True)


def normalize_bakta_type(value):
    """Return the selected Bakta database type, using full when omitted."""
    if value is None:
        return DEFAULT_BAKTA_TYPE
    if not isinstance(value, str) or value not in BAKTA_TYPES:
        raise ValueError("config.bakta.type must be either 'full' or 'light'")
    return value


def bakta_resource_subpath(selected_type):
    """Return the distinct resource-root subpath for one Bakta database type."""
    if selected_type == "full":
        return "bakta_db/db"
    if selected_type == "light":
        return "bakta_db/light/db"
    raise ValueError("Bakta database type must be either 'full' or 'light'")


def resolve_bakta_database(database_type, explicit_db, configured_root=None):
    """Resolve the selected type and explicit, environment, or default directory."""
    selected = normalize_bakta_type(database_type)
    database = resolve_resource_path(
        explicit_db,
        "BAKTA_DB",
        bakta_resource_subpath(selected),
        configured_root,
    )
    return selected, database
