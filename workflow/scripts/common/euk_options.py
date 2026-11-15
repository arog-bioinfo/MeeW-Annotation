"""Validation helpers for workflow-managed eukaryotic annotation options."""

import shlex

from .paths import expand_path

METAEUK_MANAGED_OPTIONS = {
    "--threads",
    "--max-exon-sets",
    "--exhaustive-search",
}

SUPPORTED_FUNANNOTATE2_ROUTE_MESSAGE = (
    "Funannotate2 requires one supported model route: provide params JSON "
    "(species is derived from JSON), provide pretrained plus species, or provide "
    "species for genome-specific training. Unknown-species automation using MetaEuk "
    "predictions and taxonomy-aware MMseqs2 classification is deferred and is not "
    "implemented."
)


def _managed_token_option(token, managed_options):
    """Return a managed CLI option represented by a token, including attached shorts."""
    option = token.split("=", 1)[0]
    if option in managed_options:
        return option
    for candidate in sorted(managed_options, key=len, reverse=True):
        if candidate.startswith("-") and not candidate.startswith("--"):
            if token.startswith(candidate) and token != candidate:
                return candidate
    return None


def params_json_context(path):
    """Read required species and optional taxonomy from Funannotate2 params JSON."""
    import json
    from pathlib import Path

    source = Path(expand_path(path))
    if not source.is_file():
        raise ValueError(f"Funannotate2 params JSON does not exist: {source}")
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(
            f"Unable to read Funannotate2 params JSON {source}: {error}"
        ) from error
    species = payload.get("species") if isinstance(payload, dict) else None
    if not isinstance(species, str) or not species.strip():
        raise ValueError(
            f"Funannotate2 params JSON {source} must contain a nonempty string 'species'."
        )
    taxonomy = payload.get("taxonomy")
    if taxonomy is not None and (
        not isinstance(taxonomy, dict)
        or any(
            not isinstance(key, str) or not isinstance(value, str) or not value.strip()
            for key, value in taxonomy.items()
        )
    ):
        raise ValueError(
            f"Funannotate2 params JSON {source} taxonomy must be an object of "
            "nonempty string values when provided."
        )
    return species.strip(), taxonomy


def params_json_species(path):
    """Read the required species string from a Funannotate2 params JSON."""
    return params_json_context(path)[0]


def select_funannotate2_route(configured_species, params_json, pretrained):
    """Resolve the deterministic model route and effective species."""
    species = str(configured_species).strip() if configured_species else ""
    params_json = str(params_json).strip() if params_json else ""
    pretrained = str(pretrained).strip() if pretrained else ""
    if params_json and pretrained:
        raise ValueError(
            "funannotate2.params and funannotate2.pretrained are mutually exclusive."
        )
    if params_json:
        derived_species = params_json_species(params_json)
        if species and species != derived_species:
            raise ValueError(
                "Configured Funannotate2 species disagrees with params JSON species: "
                f"configured={species!r}, params={derived_species!r}."
            )
        return "params_json", derived_species, params_json
    if pretrained:
        if not species:
            raise ValueError(
                "The Funannotate2 pretrained route requires configured species."
            )
        return "pretrained", species, pretrained
    if species:
        return "genome_specific_training", species, ""
    raise ValueError(SUPPORTED_FUNANNOTATE2_ROUTE_MESSAGE)


def validate_metaeuk_extra(extra):
    """Reject options whose values are managed by the workflow contract."""
    extra = "" if extra is None else str(extra)
    try:
        tokens = shlex.split(extra)
    except ValueError as error:
        raise ValueError(
            f"config.metaeuk.extra is not valid shell syntax: {error}"
        ) from error

    conflicts = []
    for token in tokens:
        option = _managed_token_option(token, METAEUK_MANAGED_OPTIONS)
        if option in METAEUK_MANAGED_OPTIONS:
            conflicts.append(option)
    if conflicts:
        formatted = ", ".join(repr(value) for value in sorted(set(conflicts)))
        raise ValueError(
            "config.metaeuk.extra must not set MetaEuk options managed by MeeW: "
            f"{formatted}. Configure metaeuk.max_exon_sets for that managed value; "
            "--exhaustive-search is intentionally not enabled automatically."
        )
    return extra


def validate_funannotate2_route(config, inspect_params=False):
    """Reject ambiguous/missing Funannotate2 model routes before DAG execution."""
    params = config.get("params") or ""
    pretrained = config.get("pretrained") or ""
    species = config.get("species") or ""
    if params and pretrained:
        raise ValueError(
            "config.funannotate2.params and config.funannotate2.pretrained are mutually "
            "exclusive; select exactly one reusable-model route."
        )
    if params and inspect_params:
        derived_species = params_json_species(params)
        if species and species != derived_species:
            raise ValueError(
                "Configured Funannotate2 species disagrees with params JSON species: "
                f"configured={species!r}, params={derived_species!r}."
            )
    if pretrained and not species:
        raise ValueError(
            "config.funannotate2.pretrained requires a nonempty config.funannotate2.species."
        )
    if not params and not pretrained and not species:
        raise ValueError(SUPPORTED_FUNANNOTATE2_ROUTE_MESSAGE)


def validate_funannotate2_extras(config):
    """Prevent phase extras from overriding the workflow's scientific state contract."""
    managed = {
        "extra_clean": {"-f", "--fasta", "-o", "--out", "--cpus"},
        "extra_train": {
            "-f",
            "--fasta",
            "-s",
            "--species",
            "-o",
            "--out",
            "--strain",
            "--cpus",
            "--busco-lineage",
        },
        "extra_predict": {
            "-i",
            "--input-dir",
            "-f",
            "--fasta",
            "-o",
            "--out",
            "-p",
            "--params",
            "--pretrained",
            "-s",
            "--species",
            "-st",
            "--strain",
            "-c",
            "--cpus",
        },
        "extra_annotate": {
            "-i",
            "--input-dir",
            "-f",
            "--fasta",
            "-o",
            "--out",
            "-s",
            "--species",
            "-st",
            "--strain",
            "--cpus",
            "--busco-lineage",
        },
    }
    for field, managed_options in managed.items():
        extra = config.get(field) or ""
        try:
            tokens = shlex.split(str(extra))
        except ValueError as error:
            raise ValueError(
                f"config.funannotate2.{field} is not valid shell syntax: {error}"
            ) from error
        conflicts = sorted(
            {
                _managed_token_option(token, managed_options)
                for token in tokens
                if _managed_token_option(token, managed_options)
            }
        )
        if conflicts:
            formatted = ", ".join(repr(option) for option in conflicts)
            raise ValueError(
                f"config.funannotate2.{field} must not set options managed by MeeW: "
                f"{formatted}."
            )
